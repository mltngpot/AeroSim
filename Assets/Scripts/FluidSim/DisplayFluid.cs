using System.Collections.Generic;
using UnityEngine;



namespace EulerianFluidSimulator
{
    //Display the fluid simulation data on a texture
    //Display streamlines and velocities with lines
    //Display obstacles as mesh
    public static class DisplayFluid
    {
        //Called every Update
        public static GameObject Tank;
        private static GameObject segmentTemplate;
        private static GameObject frame;

        public static void Draw(FluidScene scene)
        {
            if (Tank == null)
            {
                Tank = GameObject.Find("Tank");
                Tank.transform.localScale = new Vector3(scene.fluid.SimWidth, scene.fluid.SimHeight, scene.fluid.SimDepth);
                var collider = Tank.GetComponent<BoxCollider>();
                collider.size = new Vector3(scene.fluid.SimWidth, scene.fluid.SimHeight, scene.fluid.SimDepth);
                segmentTemplate = new GameObject("Line");
                segmentTemplate.AddComponent<LineRenderer>();
                LineRenderer lr = segmentTemplate.GetComponent<LineRenderer>();
                lr.material = new Material(Tank.GetComponent<MeshRenderer>().material);
                lr.startColor = Color.black;
                lr.endColor = Color.black;
                lr.startWidth = 0.01f;
                lr.endWidth = 0.01f;
            }

            //DrawContainer(scene);
        }

        private static void DrawContainer(FluidScene scene)
        {
            if (frame != null)
            {
                return;
            }

            FluidSim f = scene.fluid;

            Vector3 leftLowerFront = scene.SimToWorld(new Vector3(0f, 0f, 0f));
            Vector3 rightLowerFront = scene.SimToWorld(new Vector3(f.SimWidth, 0f, 0f));
            Vector3 leftUpperFront = scene.SimToWorld(new Vector3(0f, f.SimHeight, 0f));
            Vector3 rightUpperFront = scene.SimToWorld(new Vector3(f.SimWidth, f.SimHeight, 0f));
            Vector3 leftLowerBack = scene.SimToWorld(new Vector3(0f, 0f, f.SimDepth));
            Vector3 rightLowerBack = scene.SimToWorld(new Vector3(f.SimWidth, 0f, f.SimDepth));
            Vector3 leftUpperBack = scene.SimToWorld(new Vector3(0f, f.SimHeight, f.SimDepth));
            Vector3 rightUpperBack = scene.SimToWorld(new Vector3(f.SimWidth, f.SimHeight, f.SimDepth));

            frame = new GameObject("Frame");
            frame.transform.parent = Tank.transform;

            addLine(new Vector3[] { leftLowerFront, rightLowerFront }, frame);
            addLine(new Vector3[] { leftUpperFront, rightUpperFront }, frame);
            addLine(new Vector3[] { leftLowerBack, rightLowerBack }, frame);
            addLine(new Vector3[] { leftUpperBack, rightUpperBack }, frame);
            addLine(new Vector3[] { leftLowerFront, leftUpperFront }, frame);
            addLine(new Vector3[] { rightLowerFront, rightUpperFront }, frame);
            addLine(new Vector3[] { leftLowerBack, leftUpperBack }, frame);
            addLine(new Vector3[] { rightLowerBack, rightUpperBack }, frame);

        }

        private static void addLine(Vector3[] points, GameObject lineHolder)
        {
            GameObject segment = GameObject.Instantiate(segmentTemplate, lineHolder.transform);

            LineRenderer lineRenderer = segment.GetComponent<LineRenderer>();
            lineRenderer.SetPositions(points);
            lineRenderer.positionCount = points.Length;
        }



        //
        // Show the fluid simulation data on a texture
        //
        private static void UpdateTexture(FluidScene scene)
        {
            FluidSim f = scene.fluid;

            Texture3D fluidTexture = scene.fluidTexture;

            //Generate a new texture if none exists or if we have changed resolution
            if (fluidTexture == null || fluidTexture.width != f.numX || fluidTexture.height != f.numY || fluidTexture.depth != f.numZ)
            {
                fluidTexture = new(f.numX, f.numY, f.numZ, TextureFormat.RGBA32, false);

                //Dont blend the pixels
                //fluidTexture.filterMode = FilterMode.Point;

                //Blend the pixels 
                fluidTexture.filterMode = FilterMode.Bilinear;

                //Don't wrap the border with the border on the opposite side of the texture
                fluidTexture.wrapMode = TextureWrapMode.Clamp;

                scene.fluidMaterial.SetTexture("texture", fluidTexture);

                scene.fluidTexture = fluidTexture;
            }


            //The texture colors
            Color32[] textureColors = new Color32[f.numX * f.numY * f.numZ];

            //Find min and max pressure
            MinMax minMaxP = f.GetMinMaxPressure();

            //Find the colors
            //This was an array in the source, but we can treat the Vector4 as an array to make the code match
            //Vector4 color = new (255, 255, 255, 255);

            for (int i = 0; i < 1; i++)
            {
                for (int j = 0; j < f.numY; j++)
                {
                    for (int k = 0; k < f.numZ; k++)
                    {
                        //This was an array in the source, but we can treat the Vector4 as an array to make the code match
                        //Moved to here from before the loop so it resets every time so we can display the walls if we deactivate both pressure and smoke
                        Vector4 color = new(255, 255, 255, 255);

                        if (scene.showPressure)
                        {
                            float p = f.p[f.To1D(i, j, k)];

                            //Blue means low pressure and red is high pressure
                            color = UsefulMethods.GetSciColor(p, minMaxP.min, minMaxP.max);

                            //Color the smoke according to the scientific color scheme 
                            //Everything that's not smoke becomes black
                            //Everything that's smoke shows the pressure field
                            if (scene.showSmoke)
                            {
                                //How much smoke in this cell?
                                float smoke = f.m[f.To1D(i, j, k)];

                                //smoke = 0 means max smoke, so will be 0 if no smoke in the cell (smoke = 1)
                                color[0] = Mathf.Max(0f, color[0] - 255 * smoke);
                                color[1] = Mathf.Max(0f, color[1] - 255 * smoke);
                                color[2] = Mathf.Max(0f, color[2] - 255 * smoke);
                            }
                        }
                        else if (scene.showSmoke)
                        {
                            //How much smoke in this cell?
                            float smoke = f.m[f.To1D(i, j, k)];

                            //smoke = 0 means max smoke, and 255 * 0 = 0 -> black 
                            color[0] = 255 * smoke;
                            color[1] = 255 * smoke;
                            color[2] = 255 * smoke;

                        }
                        //If both pressure and smoke are deactivated, then display obstacles as black, the rest as white
                        //There was a bug in the source code where everything turned back, but "f.s[f.To1D(i, j)] == 0f" should mean only walls should be black
                        else if (f.s[f.To1D(i, j, k)] == 0f)
                        {
                            color[0] = 0;
                            color[1] = 0;
                            color[2] = 0;
                        }

                        //Add the color to this pixel
                        //Color32 is 0-255
                        Color32 pixelColor = new((byte)color[0], (byte)color[1], (byte)color[2], (byte)color[3]);

                        textureColors[f.To1D(i, j, k)] = pixelColor;
                    }
                }
            }

            //Add all colors to the texture
            fluidTexture.SetPixels32(textureColors);

            //Copies changes you've made in a CPU texture to the GPU
            fluidTexture.Apply(false);
        }



        //
        // Show the u and v velocities at each cell by drawing lines
        //
        private static void ShowVelocities(FluidScene scene)
        {
            FluidSim f = scene.fluid;

            //Cell width
            float h = f.h;

            //The length of the lines which will be scaled by the velocity in simulation space
            float scale = 0.02f;



            GameObject velocities = new GameObject("Velocities");
            velocities.transform.parent = Tank.transform;

            float maxDimension = Mathf.Max(f.numX, f.numY, f.numZ);
            float[] floats = new float[(int)Mathf.Ceil(maxDimension) + 1];
            for (int i = 0; i < floats.Length; i++)
            {
                floats[i] = i * h;
            }

            for (int i = 0; i < f.numX; i++)
            {
                Debug.Log(i);
                for (int j = 0; j < f.numY; j++)
                {
                    for (int k = 0; k < f.numZ; k++)
                    {
                        float u = f.u[f.To1D(i, j, k)];
                        float v = f.v[f.To1D(i, j, k)];
                        float w = f.w[f.To1D(i, j, k)];


                        float x1 = floats[i] + u * scale;
                        float y1 = floats[j] + v * scale;
                        float z1 = floats[k] + w * scale;

                        Vector3 start = new Vector3(floats[i], floats[j], floats[k]);
                        Vector3 end = new Vector3(x1, y1, z1);

                        addLine(new Vector3[] { start, end}, velocities);
                    }
                }
            }
        }

        //
        // Show streamlines that follows the velocity to easier visualize how the fluid flows
        //
        private static void ShowStreamlines(FluidScene scene)
        {
            FluidSim f = scene.fluid;

            //How many segments per streamline?
            int numSegs = 15;

            List<Vector3> streamlineCoordinates = new();

            GameObject streamlines = new GameObject("Streamlines");
            streamlines.transform.parent = Tank.transform;

            //Dont display a streamline from each cell because it makes it difficult to see, so every 5 cell
            for (int i = 1; i < f.numX - 1; i += 5)
            {
                for (int j = 1; j < f.numY - 1; j += 5)
                {
                    for (int k = 1; k < f.numZ - 1; k += 5)
                    {
                        //Reset
                        streamlineCoordinates.Clear();

                        //Center of the cell in simulation space
                        Vector3 samplePosition = new(i * f.h, j * f.h, k * f.h);

                        //Simulation space to global
                        Vector3 startPos = scene.SimToWorld(samplePosition);

                        streamlineCoordinates.Add(startPos);

                        //Build the line
                        for (int n = 0; n < numSegs; n++)
                        {
                            //The velocity at the current coordinate
                            float u = f.SampleField(samplePosition, FluidSim.SampleArray.uField);
                            float v = f.SampleField(samplePosition, FluidSim.SampleArray.vField);
                            float w = f.SampleField(samplePosition, FluidSim.SampleArray.wField);

                            //Move a small step in the direction of the velocity
                            samplePosition.x += u * 0.01f;
                            samplePosition.y += v * 0.01f;
                            samplePosition.z += w * 0.01f;

                            //Stop the line if we are outside of the simulation area
                            //The guy in the video is only checking x > f.GetWidth() for some reason...
                            if (samplePosition.x > f.SimWidth || samplePosition.x < 0f ||
                                samplePosition.y > f.SimHeight || samplePosition.y < 0f ||
                                samplePosition.z > f.SimDepth || samplePosition.z < 0f)
                            {
                                break;
                            }

                            //Add the next coordinate of the streamline

                            streamlineCoordinates.Add(samplePosition);
                        }
                        try
                        {
                            //Display the line
                            addLine(streamlineCoordinates.ToArray(), streamlines);
                        }
                        catch (System.Exception e)
                        {
                            Debug.Log(e);
                        }
                    }
                }
            }
        }
    }
}
