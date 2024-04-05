using UnityEngine;

//General class to interpolate between points in a grid
//These points can be in a staggered grid where they are defined on the horizontal and vertical lines
public static class GridInterpolation
{
    //Derivation of how to find the linear interpolation of P by using A, B, C, D and their respective coordinates
    //This is a square with side length h (I did my best)
    // C------D
    // |      |
    // |___P  |
    // |   |  |
    // A------B
    // The points have the coordinates
    // A: (xA, yA)
    // B: (xB, yB)
    // C: (xC, yC)
    // D: (xD, yD)
    // P: (xP, yP)
    //
    //We need to do 3 linear interpolations to find P:
    // P_AB = (1 - tx) * A + tx * B
    // P_CD = (1 - tx) * C + tx * D
    // P = (1 - ty) * P_AB + ty * P_CD
    //
    //Insert P_AB and P_CD into P and we get: 
    // P = (1 - ty) * [(1 - tx) * A + tx * B] + ty * [(1 - tx) * C + tx * D] 
    // P = (1 - tx) * (1 - ty) * A + tx * (1 - ty) * B + (1 - tx) * ty * C + tx * ty * D
    //
    //t is a parameter in the range [0, 1]. If tx = 0 we get A or if tx = 1 we get B in the P_AB case
    //The parameter can be written as:
    // tx = (xp - xA) / (xB - xA) = (xP - xA) / h = deltaX / h
    // ty = (yp - yA) / (yB - yA) = (yP - yA) / h = deltaY / h
    //
    //Define:
    // sx = 1 - tx
    // sy = 1 - ty
    //
    //And we get the following:
    // P = sx * sy * A + tx * sy * B + sx * ty * C + tx * ty * D
    //
    //Simplify the weights:
    // wA = sx * sy
    // wB = tx * sy
    // wC = sx * ty
    // wD = tx * ty
    //
    //Note that: wA + wB + wC + wD = 1
    //
    //The final iterpolation:
    // P = wA * A + wB * B + wC * C + wD * D 
    //
    //In simple code (which is slightly slower than the above because we do some calculations multiple times but easy to understand):
    //float tx = Mathf.InverseLerp(xA, xB, xP);
    //float ty = Mathf.InverseLerp(yA, yB, yP);
    //float P_AB = Mathf.Lerp(A, B, tx); 
    //float P_CD = Mathf.Lerp(C, D, tx);
    //float P = Mathf.Lerp(P_AB, P_CD, ty);
    public static GridWeights GetWeights(this GridConstants gridData, Vector3 P, Vector3 A)
    {
        float deltaX = P.x - A.x;
        float deltaY = P.y - A.y;
        float deltaZ = P.z - A.z;

        float tx = deltaX * gridData.one_over_h;
        float ty = deltaY * gridData.one_over_h;
        float tz = deltaZ * gridData.one_over_h;

        float sx = 1 - tx;
        float sy = 1 - ty;
        float sz = 1 - tz;

        GridWeights weights = new GridWeights
        {
            A = sx * sy * sz,
            B = tx * sy * sz,
            C = sx * ty * sz,
            D = tx * ty * sz,
            E = sx * sy * tz,
            F = tx * sy * tz,
            G = sx * ty * tz,
            H = tx * ty * tz
        };

        return weights;
    }



    //3x3 Staggered grid (c is not staggered)
    //0,0 is bottom left
    //Notice there's no u on the vertical line to the right of the last cell and no v at the top
    //+-----+-----+-----+
    //|     |     |     |
    //u  c  u  c  u  c  |
    //|     |     |     |
    //+--v--+--v--+--v--+
    //|     |     |     |
    //u  c  u  c  u  c  |
    //|     |     |     |
    //+--v--+--v--+--v--+
    //|     |     |     |
    //u  c  u  c  u  c  |
    //|     |     |     |
    //+--v--+--v--+--v--+

    //The values we want to interpolate can be on:
    //- The center of the cell
    //- In the middle of the vertical lines (staggered grid): u
    //- In the middle of the horizontal lines (staggered grid): v
    public enum Grid
    {
        u,
        v,
        w,
        center
    }



    //Clamp the interpolation point P so we know we can interpolate from 4 grid points
    public static Vector3 ClampInterpolationPoint(this GridConstants gridData, Vector3 P, Grid sampleField)
    {
        Vector3 minOffset;
        Vector3 maxOffset;

        //Which grid to we want to interpolate from? 
        switch (sampleField)
        {
            case Grid.u:
                minOffset = new Vector3(0f, gridData.half_h, gridData.half_h);
                maxOffset = new Vector3(gridData.h, gridData.half_h, gridData.half_h);
                break;

            case Grid.v:
                minOffset = new Vector3(gridData.half_h, 0f, gridData.half_h);
                maxOffset = new Vector3(gridData.half_h, gridData.h, gridData.half_h);
                break;

            case Grid.w:
                minOffset = new Vector3(gridData.half_h, gridData.half_h, 0f);
                maxOffset = new Vector3(gridData.half_h, gridData.half_h, gridData.h);
                break;

            case Grid.center:
            default:
                minOffset = new Vector3(gridData.half_h, gridData.half_h, gridData.half_h);
                maxOffset = new Vector3(gridData.half_h, gridData.half_h, gridData.half_h);
                break;
        }

        //In the center case we have to push it in 0.5*h in all directions
        //In the u case:
        //- We have to push it in 0 from the left side because the point can be on the vertical line
        //- h from the right side because there's no u on the right side of the last cell
        //In the v case we have to push it in 0.5 from both left and right side in x direction
        Vector3 clamped = new Vector3(
                       Mathf.Min(P.x, gridData.gridSize.x - maxOffset.x),
                       Mathf.Min(P.y, gridData.gridSize.y - maxOffset.y),
                       Mathf.Min(P.z, gridData.gridSize.z - maxOffset.z));

        //Clamp to be inside of the grid on the minimum side
        clamped = new Vector3(
                       Mathf.Max(clamped.x, minOffset.x),
                       Mathf.Max(clamped.y, minOffset.y),
                       Mathf.Max(clamped.z, minOffset.z));

        return clamped;
    }



    //Get grid indices of A which is defined to always be to the left and below the sample point 
    //We assume the interpolation point has been clamped
    //We only need the indices of A, the other ones are just xA_index + 1 and yA_index + 1 
    public static Vector3Int GetAIndices(this GridConstants gridData, Vector3 P, Grid sampleField)
    {
        Vector3 d = GetGridOffsets(sampleField, gridData.half_h);

        //To go from coordinate to cell we generally do: FloorToInt(pos / cellSize)
        //But we have to compensate because of the staggered grid to always get the index to the "left" of the point (in x direction) 
        //If we want to sample the center of the cell then dx will be half the cell width to push the sample point into the correct cell to always get the index of A which is to the left of the sample point 
        int cellIndexX = Mathf.FloorToInt((P.x - d.x) * gridData.one_over_h);
        int cellIndexY = Mathf.FloorToInt((P.y - d.y) * gridData.one_over_h);
        int cellIndexZ = Mathf.FloorToInt((P.z - d.z) * gridData.one_over_h);

        //Clamp
        //The input position is assumed to be within the grid
        //To interpolate between 4 values we have to make sure we are 1 extra cell in from the maximum side
        //So on the maximum side we will interpolate between gridData.numX - 2 and gridData.numX - 1 in x direction 

        Vector3Int indicies = Vector3Int.zero;
        indicies.x = Mathf.Max(0, Mathf.Min(cellIndexX, gridData.arrayLengths.x - 2));
        indicies.y = Mathf.Max(0, Mathf.Min(cellIndexY, gridData.arrayLengths.y - 2));
        indicies.z = Mathf.Max(0, Mathf.Min(cellIndexZ, gridData.arrayLengths.z - 2));

        return indicies;
    }



    //Get grid coordinates of A in local space
    //We only need the coordinates of A, the other ones are just xA + cellWidth and yA + cellWidth 
    public static Vector3 GetACoordinates(this GridConstants gridData, Grid sampleField, Vector3Int indicies)
    {
        Vector3 d = GetGridOffsets(sampleField, gridData.half_h);
        Vector3 A = Vector3.zero;

        A.x = indicies.x * gridData.h + d.x;
        A.y = indicies.y * gridData.h + d.y;
        A.z = indicies.z * gridData.h + d.z;

        return A;
    }



    //Get parameters so we can standardize the code depending on which grid data we want to sample
    //+-----+-----+-----+
    //|     |     |     |
    //u  c  u  c  u  c  |
    //|     |     |     |
    //+--v--+--v--+--v--+
    //|     |     |     |
    //u  c  u  c  u  c  |
    //|     |     |     |
    //+--v--+--v--+--v--+
    //|     |     |     |
    //u  c  u  c  u  c  |
    //|     |     |     |
    //+--v--+--v--+--v--+
    //If we want to sample the center of the cell then dx will be half the cell width to push the sample point into the correct cell to always get the index of A which is to the left of the sample point in x direction
    private static Vector3 GetGridOffsets(Grid sampleField, float half_h)
    {
        Vector3 d = Vector3.zero;

        //Which grid to we want to interpolate from? 
        switch (sampleField)
        {
            case Grid.u: 
                d.y = half_h; 
                break;

            case Grid.v: 
                d.x = half_h; 
                break;

                case Grid.w:
                    d.z = half_h;
                break;

            case Grid.center: 
                d.x = half_h; 
                d.y = half_h; 
                d.z = half_h;
                break;
        }

        return d;
    }
}
