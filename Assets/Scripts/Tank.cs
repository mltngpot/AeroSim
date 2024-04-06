using System;
using System.Reflection;
using System.Security.Cryptography;
using Unity.VisualScripting;
using UnityEngine;
using UnityEngine.UIElements;

public class Tank : MonoBehaviour
{
    [Header("Simulation Settings")]
    public float overRelaxation = 1.9f;
    public uint resolution = 10;
    public uint incompessibilityIterrations = 100;
    public float Density = 1000f;
    public Vector3 externalForce = new Vector3(0, 0, 0);
    public Vector3 wind = new Vector3(0, 0, 0);
    public float deltaTime = 1 / 60f;

    [Header("Tank Settings")]
    public GameObject cellPrefab;
    public bool showVelocity = true;
    public bool showStreamLine = false;
    public bool showPressure = false;
    public bool showSmoke = false;

    [Header("Shader Settings")]
    public ComputeShader shader;
    ComputeBuffer buffer;
    private int integrationKernel;
    private int incompressibilityKernel;
    private int extrapolationKernel;
    private int preAdvectionUpdateKernel;
    private int advectionKernel;
    private int postAdvectionUpdateKernel;

    private Cell[] cells;
    private CellData[] cellData;
    private Vector3Int tankSize;
    private float h;
    private float half_h;
    private float inverted_h;
    private int numCells;
    private int numCellData;
    private bool unPaused = false;
    private float computedPressure;
    private int shaderGroups;

    public void Pause()
    {
        unPaused = !unPaused;
    }

    private void Start()
    {
        BuildTankConstants();
        AssembleCells();
        InitializeShader();
    }

    private void BuildTankConstants()
    {
        tankSize = new Vector3Int((int)(transform.localScale.x * resolution), (int)(transform.localScale.y * resolution), (int)(transform.localScale.z * resolution));
        h = 1 / resolution;
        half_h = h / 2;
        inverted_h = 1 / h;
        computedPressure = Density * h / deltaTime;
        numCells = (tankSize.x) * (tankSize.y)  * (tankSize.z);
        cells = new Cell[numCells];
        numCellData = (tankSize.x + 2) * (tankSize.y + 2) * (tankSize.z + 2);
        cellData = new CellData[numCellData];

        for (int index = 0; index < numCellData; index++)
        {
            CellData cellData = new CellData();
            Vector3Int position = CellDataToPosition(index);
            if (position.x == 0 || position.y == 0 || position.z == 0 || position.y == tankSize.y + 1)
            {
                cellData.permeability = 0f;
            }
            else
            {
                cellData.permeability = 1f;
            }
            if (position.x == 0)
            {
                cellData.velocity.x = wind.x;
            }
            if (position.z == 0)
            {
                cellData.velocity.z = wind.z;
            }

            this.cellData[index] = cellData;
        }
    }

    private Vector3Int CellToPosition(int index)
    {
        int x = index / (tankSize.y * tankSize.z);
        int y = (index - x * tankSize.y * tankSize.z) / tankSize.z;
        int z = index - x * tankSize.y * tankSize.z - y * tankSize.z;
        return new Vector3Int(x, y, z);
    }

    private Vector3Int CellDataToPosition(int index)
    {
        int x = index / ((tankSize.y + 2) * (tankSize.z + 2));
        int y = (index - x * (tankSize.y + 2) * (tankSize.z + 2)) / (tankSize.z + 2);
        int z = index - x * (tankSize.y + 2) * (tankSize.z + 2) - y * (tankSize.z + 2);
        return new Vector3Int(x, y, z);
    }

    private int CellTo1D(int x, int y, int z)
    {
        return ((x * tankSize.y) + y) * tankSize.z + z;
    }

    private int CellDataTo1D(int x, int y, int z)
    {
        return ((x * (tankSize.y + 2)) + y) * (tankSize.z + 2) + z;
    }

    private int CellIndexToCellDataIndex(int cellIndex)
    {
        var position = CellToPosition(cellIndex);
        return CellDataTo1D(position.x + 1, position.y + 1, position.z + 1);
    }

    private void AssembleCells()
    {
        buildTankGameObjects();
        Cell[] cells = gameObject.GetComponentsInChildren<Cell>();
        for (int c = 0; c < cells.Length; c++)
        {
            Cell cell = cells[c];
            Vector3 position = cell.transform.localPosition +
                               cell.transform.parent.localPosition +
                               cell.transform.parent.parent.localPosition;
            int index = CellTo1D((int)position.x, (int)position.y, (int)position.z);
            int dataIndex = CellIndexToCellDataIndex(index);
            cell.data = cellData[dataIndex];
            cell.setId(position);
            this.cells[index] = cell;
        }
    }

    private void buildTankGameObjects()
    {
        GameObject row = new GameObject("Row");
        for (int i = 0; i < tankSize.x; i++)
        {
            GameObject cell = Instantiate(cellPrefab, row.transform);
            cell.name = $"Cell_{i}";
            cell.transform.localPosition = new Vector3(i, 0, 0);
            cell.transform.parent = row.transform;
        }
        GameObject plane = new GameObject("Plane");
        for (int j = 0; j < tankSize.y; j++)
        {
            GameObject rowInstance = Instantiate(row, plane.transform);
            rowInstance.name = $"Row_{j}";
            rowInstance.transform.localPosition = new Vector3(0, j, 0);
            rowInstance.transform.parent = plane.transform;
        }
        for (int k = 0; k < tankSize.z; k++)
        {
            GameObject planeInstance = Instantiate(plane, gameObject.transform);
            planeInstance.name = $"Plane_{k}";
            planeInstance.transform.localPosition = new Vector3(0, 0, k);
            planeInstance.transform.parent = gameObject.transform;
        }
        gameObject.transform.localScale = new(h, h, h);
        Destroy(row);
        Destroy(plane);
    }
    
    private void FixedUpdate()
    {
        if (unPaused)
        {
            unPaused = false;
            Simulate();
            UpdateCells();
        }
    }

    private void UpdateCells()
    {
        buffer.GetData(cellData);
        for (int i = 0; i < tankSize.x; i++)
        {
            for (int j = 0; j < tankSize.y; j++)
            {
                for (int k = 0; k < tankSize.z; k++)
                {
                    int index = CellTo1D(i, j, k);
                    try
                    {
                        cells[index].UpdateCell(cellData[CellIndexToCellDataIndex(index)]);
                        cells[index].showVelocity = this.showVelocity;
                        cells[index].showStreamLine = this.showStreamLine;
                    }
                    catch (Exception e)
                    {
                        Debug.Log($"({i}, {j}, {k}) {index}");
                    }
                }
            }
        }
    }

    private void InitializeShader()
    {
        int cellSize = sizeof(float) * 10;
        buffer = new ComputeBuffer(cellData.Length, cellSize);

        integrationKernel = shader.FindKernel("Integration");
        incompressibilityKernel = shader.FindKernel("SolveIncompressability");
        extrapolationKernel = shader.FindKernel("Extrapolation");
        preAdvectionUpdateKernel = shader.FindKernel("PreAdvectionUpdate");
        advectionKernel = shader.FindKernel("Advection");
        postAdvectionUpdateKernel = shader.FindKernel("PostAdvectionUpdate");

        shader.SetBuffer(integrationKernel, "cellData", buffer);
        shader.SetBuffer(incompressibilityKernel, "cellData", buffer);
        shader.SetBuffer(extrapolationKernel, "cellData", buffer);
        shader.SetBuffer(preAdvectionUpdateKernel, "cellData", buffer);
        shader.SetBuffer(advectionKernel, "cellData", buffer);
        shader.SetBuffer(postAdvectionUpdateKernel, "cellData", buffer);

        shader.SetInt("tankSizeX", tankSize.x);
        shader.SetInt("tankSizeY", tankSize.y);
        shader.SetInt("tankSizeZ", tankSize.z);
        shader.SetInt("XStep", tankSize.y * tankSize.z);
        shader.SetInt("YStep", tankSize.z);
        shader.SetInt("ZStep", 1);
        shader.SetFloat("overRelaxation", overRelaxation);
        shader.SetFloat("computedPressure", computedPressure);
        shader.SetFloat("tankScale", h);
        shader.SetFloat("halfTankScale", half_h);
        shader.SetFloat("invertedTankScale", inverted_h);
        shader.SetFloat("deltaTime", Time.fixedDeltaTime);

        shaderGroups = buffer.count / 32;
    }


    private void Simulate()
    {
        Integrate();
        SolveForIncompressability();
        buffer.GetData(cellData);
        Extrapolate();
        buffer.GetData(cellData);
        Advection();
    }

    private void Integrate()
    {
        //shader.Dispatch(integrationKernel, 64, 1, 1);
        //doing nothing for now
    }

    private void SolveForIncompressability()
    {
        for (int i = 0; i < incompessibilityIterrations; i++)
        {
            shader.Dispatch(incompressibilityKernel, shaderGroups, 1, 1);
        }
    }

    private void Extrapolate()
    {
        shader.Dispatch(extrapolationKernel, shaderGroups, 1, 1);
    }

    private void Advection()
    {
        shader.Dispatch(preAdvectionUpdateKernel, shaderGroups, 1, 1);
        shader.Dispatch(advectionKernel, shaderGroups, 1, 1);
        shader.Dispatch(postAdvectionUpdateKernel, shaderGroups, 1, 1);
    }
}
