using System.Collections;
using System.Collections.Generic;
using UnityEngine;
using EulerianFluidSimulator;
using System;
using System.Linq.Expressions;
using UnityEditor;

//2D Fluid Simulator
//Based on: "How to write an Eulerian Fluid Simulator with 200 lines of code" https://matthias-research.github.io/pages/tenMinutePhysics/
//You can pause the simulation with P and when paused you can move the simulation in steps with M
//Eulerian means we simulate the fluid in a grid - not by using particles (Lagrangian)
//Can simulate both liquids and gas. But we will always use water density because we only measure the pressure distribution in the water tank. Density is only affecting the pressure calculations - not the velocity field, so it doesn't matter
//Assume incompressible fluid with zero viscosity (inviscid) which are good approximations for water and gas
//Known bugs from the original code (some have been fixed):
// - Integrate() is ignoring the last column in x direction
// - The wind tunnel in-velocities are set to zero if we move the obstacle across the first columns because how the move obstacle method works
// - In Extrapolate() we should detect if it's an obstacle, so now we copy velocities into obstacles. The result is the same but can be confusing
// - In "paint" the smoke is not reaching within 2 cells on the right and top side
// - When we deactivate both pressure and smoke, we want to display the obstacles/walls as black and the fluid as white. In the source code there was a bug where everything turned black even though (at least I think so based on the source code) that only the walls should be black. This was fixed in my code
// - The SampeField method has a bug in it where if we try to sample from the bottom row we interpolate from the from values. This was fixed in my code

public class FluidSimController : MonoBehaviour
{
    //Public
    public Material fluidMaterial;


    //Private
    private FluidScene scene;

    private FluidUI fluidUI;



    private void Start()
    {
        scene = new FluidScene(fluidMaterial);

        fluidUI = new FluidUI(this);

        SetupScene();
    }

    private void FixedUpdate()
    {
        if (Input.anyKey)
        {
            Transform cameraTransform = Camera.main.transform;
            Vector3 moveVector = new Vector3(0, 0, 0);
            Vector3 rotationAngle = new Vector3(0, 0, 0);

            if (Input.GetKey(KeyCode.W))
            {
                moveVector += cameraTransform.forward;
            }
            if (Input.GetKey(KeyCode.S))
            {
                moveVector -= cameraTransform.forward;
            }
            if (Input.GetKey(KeyCode.Q))
            {
                moveVector -= cameraTransform.right;
            }
            if (Input.GetKey(KeyCode.E))
            {
                moveVector += cameraTransform.right;
            }
            if (Input.GetKey(KeyCode.R))
            {
                moveVector += cameraTransform.up;
            }
            if (Input.GetKey(KeyCode.F))
            {
                moveVector -= cameraTransform.up;
            }
            if (Input.GetKey(KeyCode.A))
            {
                rotationAngle.y -= 10;
            }
            if (Input.GetKey(KeyCode.D))
            {
                rotationAngle.y += 10;
            }
            if (Input.GetKey(KeyCode.T))
            {
                rotationAngle.x -= 10;
            }
            if (Input.GetKey(KeyCode.G))
            {
                rotationAngle.x += 10;
            }
            Vector3 currentPos = cameraTransform.position;
            Camera.main.transform.position = new Vector3(currentPos.x + moveVector.x, currentPos.y + moveVector.y, currentPos.z + moveVector.z);
            Camera.main.transform.Rotate(rotationAngle.x, rotationAngle.y, rotationAngle.z);
            if (Input.GetKeyDown(KeyCode.Space))
            {
                Debug.Log("Simulating");
                Simulate();
                Debug.Log("Simulated");
                DisplayFluid.Draw(scene);
                Debug.Log("Drawn");
            }
        }
    }



    private void OnGUI()
    {
        fluidUI.MyOnGUI(scene);
    }



    //Simulate the fluid
    //Needs to be accessed from the UI so we can simulate step by step by pressing a key
    public void Simulate()
    {
        if (!scene.isPaused)
        {
            scene.fluid.Simulate(scene.dt, scene.gravity, scene.numIters, scene.overRelaxation);

            scene.frameNr++;
        }
    }



    //
    // Init a specific fluid simulation
    //
    public void SetupScene()
    {
        scene.overRelaxation = 1.9f;

        scene.SetTimeStep(1f / 60f);
        scene.numIters = 40;

        int res = 10;

        //The size of a cell
        float h = res;

        //How many cells do we have
        //y is up
        int numY = res;
        //The plane we use here is twice as wide as high
        int numX = numY;
        int numZ = numY;

        //Density of the fluid (water)
        float density = 1000f;

        //The dimensions of the simulation
        Vector3Int dimensions = new Vector3Int(numX, numY, numZ);

        //Create a new fluid simulator
        FluidSim f = new FluidSim(density, dimensions, h);
        scene.fluid = f;



        //The size of the plane we run the simulation on so we can convert from world space to simulation space
        scene.simTankWidth = res;
        scene.simTankHeight = res;
        scene.simTankDepth = res;

        scene.tankOffset = new(scene.simTankWidth * 0.5f, scene.simTankHeight * 0.5f, scene.simTankDepth * 0.5f);
        scene.tankScale = new(scene.fluid.SimWidth / scene.simTankWidth, scene.fluid.SimHeight / scene.simTankHeight, scene.fluid.SimDepth / scene.simTankDepth);

        //Init the different simulations
        //Wind velocity
        float inVel = 2f;

        //Set which cells are fluid or wall (default is wall = 0)
        //Also add the velocity
        for (int i = 0; i < f.numX; i++)
        {
            for (int j = 0; j < f.numY; j++)
            {
                for (int k = 0; k < f.numZ; k++)
                {
                    //1 means fluid
                    float s = 1f;

                    //Left wall, bottom wall, top wall
                    if (i == 0 || j == 0 || j == f.numY - 1 || k == 0 || k == f.numZ - 1)
                    //No right wall because we need outflow from the wind tunnel
                    //if (i == 0 || j == 0 || j == f.numY - 1 || i == f.numX - 1)
                    //Left wall
                    //if (i == 0)
                    {
                        //0 means solid
                        s = 0f;
                    }

                    f.s[f.To1D(i, j, k)] = s;

                    //Add right velocity to the fluid in the second column
                    //We now have a velocity from the wall
                    //Velocities from walls can't be modified in the simulation
                    if (i == 1)
                    {
                        f.u[f.To1D(i, j, k)] = inVel;
                    }
                }
            }
        }

        //Add smoke
        float pipeH = 0.1f * f.numY;

        //In the middle of the simulation height in y direction
        int minJ = Mathf.FloorToInt(0.5f * f.numY - 0.5f * pipeH);
        int maxJ = Mathf.FloorToInt(0.5f * f.numY + 0.5f * pipeH);
        int minK = Mathf.FloorToInt(0.5f * f.numZ - 0.5f * pipeH);
        int maxK = Mathf.FloorToInt(0.5f * f.numZ + 0.5f * pipeH);

        for (int j = minJ; j < maxJ; j++)
        {
            for (int k = minK; k < maxK; k++)
            {
                //Add the smoke in the center of the first column
                int i = 0;

                //0 means max smoke
                f.m[f.To1D(i, j, k)] = 0f;
            }
        }

        scene.gravity = 0f; //Adding gravity will break the smoke
        scene.showPressure = false;
        scene.showSmoke = true;
        scene.showStreamlines = false;
        scene.showVelocities = false;

        scene.SetTimeStep(1f / 120f);
        scene.numIters = 100;
        scene.showPressure = true;
    }
}
