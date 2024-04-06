using System.Runtime.InteropServices;
using UnityEngine;

[System.Serializable]
[StructLayout(LayoutKind.Sequential, Size = 40)]
public struct CellData
{
    public Vector3 velocity; // 12
    public Vector3 newVelocity; // 24
    public float permeability; // 28
    public float pressure; // 32
    public float density; // 36
    public float newDensity; // 40
}