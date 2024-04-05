using UnityEngine;



namespace EulerianFluidSimulator
{
    //User interactions with the fluid
    //Buttons and checkboxes
    //Position the obstacle with the mouse
    //Pause simulation (P) and step forward the simulation (M)
    //Sample cells with mouse position
    public class FluidUI
    {
        private readonly FluidSimController controller;

        public FluidUI(FluidSimController controller)
        {
            this.controller = controller;
        }

        //Buttons, checkboxes, show min/max pressure
        public void MyOnGUI(FluidScene scene)
        {
            GUILayout.BeginHorizontal("box");

            int fontSize = 20;

            RectOffset offset = new(5, 5, 5, 5);


            //Buttons
            GUIStyle buttonStyle = new(GUI.skin.button)
            {
                //buttonStyle.fontSize = 0; //To reset because fontSize is cached after you set it once 

                fontSize = fontSize,
                margin = offset
            };

            //Checkboxes
            GUIStyle toggleStyle = GUI.skin.GetStyle("Toggle");

            toggleStyle.fontSize = fontSize;
            toggleStyle.margin = offset;

            scene.showStreamlines = GUILayout.Toggle(scene.showStreamlines, "Streamlines", toggleStyle);

            scene.showVelocities = GUILayout.Toggle(scene.showVelocities, "Velocities");

            scene.showPressure = GUILayout.Toggle(scene.showPressure, "Pressure");

            scene.showSmoke = GUILayout.Toggle(scene.showSmoke, "Smoke");

            scene.useOverRelaxation = GUILayout.Toggle(scene.useOverRelaxation, "Overrelax");

            scene.overRelaxation = scene.useOverRelaxation ? 1.9f : 1.0f;

            GUILayout.EndHorizontal();


            if (scene.fluid == null)
            {
                return;
            }

            //Find min and max pressure
            MinMax minMaxP = scene.fluid.GetMinMaxPressure();

            int intMinP = Mathf.RoundToInt(minMaxP.min);
            int intMaxP = Mathf.RoundToInt(minMaxP.max);

            string pressureText = $"Pressure: {intMinP}, {intMaxP} N/m";

            GUIStyle textStyle = GUI.skin.GetStyle("Label");

            textStyle.fontSize = fontSize;
            textStyle.margin = offset;

            GUILayout.Label(pressureText, textStyle);
        }
    }
}