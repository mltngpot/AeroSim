using UnityEngine;

public record GridConstants
{
    //Cellsize
    public float h;
    public float half_h;
    public float one_over_h;

    //Grid size
    public Vector3Int gridSize;
    public Vector3Int arrayLengths;


    public GridConstants(float h, int numX, int numY, int numZ) : this(h, new Vector3Int(numX, numY, numZ)) { }

    public GridConstants(float h, Vector3Int gridSize)
    {
        this.h = h;
        half_h = h * 0.5f;
        one_over_h = 1.0f / h;

        this.gridSize = gridSize;
        arrayLengths = gridSize;
    }
}
