using System;
using UnityEngine;

public class Movement : MonoBehaviour
{
    //float rotationX = 0f;
    //float rotationY = 0f;

    public float sensitivity = 15f;

    void Update()
    {
        //rotationY += Input.GetAxis("Mouse X") * sensitivity;
        //rotationX += Input.GetAxis("Mouse Y") * -1 * sensitivity;
        //transform.localEulerAngles = new Vector3(rotationX, rotationY, 0);
    }

    private void FixedUpdate()
    {
        Move();
    }

    private void Move()
    {
        Transform cameraTransform = Camera.main.transform;
        Vector3 moveVector = new Vector3(0, 0, 0);

        if (Input.GetKey(KeyCode.W))
        {
            moveVector += cameraTransform.forward;
        }
        if (Input.GetKey(KeyCode.S))
        {
            moveVector -= cameraTransform.forward;
        }
        if (Input.GetKey(KeyCode.A))
        {
            moveVector -= cameraTransform.right;
        }
        if (Input.GetKey(KeyCode.D))
        {
            moveVector += cameraTransform.right;
        }
        if (Input.GetKey(KeyCode.Q))
        {
            moveVector += cameraTransform.up;
        }
        if (Input.GetKey(KeyCode.E))
        {
            moveVector -= cameraTransform.up;
        }

        Vector3 currentPos = cameraTransform.position;
        moveVector.Scale(new Vector3(.1f, .1f, .1f));
        Camera.main.transform.position = new Vector3(currentPos.x + moveVector.x, currentPos.y + moveVector.y, currentPos.z + moveVector.z);
    }
}
