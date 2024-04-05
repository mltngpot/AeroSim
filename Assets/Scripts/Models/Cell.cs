using UnityEngine;

public class Cell : MonoBehaviour
{
    public Material lineMaterial;
    public CellData data;
    public Vector3 id { get; private set; }
    public GameObject StreamLine;
    public GameObject Velocity;
    public bool showVelocity = true;
    public bool showStreamLine = true;

    public Cell(Vector3 id)
    {
        this.id = id;
    }

    public void setId(Vector3 id)
    {
        this.id = id;
    }
    // Start is called before the first frame update
    void Start()
    {

    }

    // Update is called once per frame
    void Update()
    {

    }

    private void FixedUpdate()
    {
        UpdateVelocity();
        UpdateStreamLine();
    }

    private void UpdateStreamLine()
    {
        if (showStreamLine)
        {
            StreamLine.GetComponent<LineRenderer>().SetPosition(0, Vector3.zero);
            StreamLine.GetComponent<LineRenderer>().SetPosition(1, data.velocity);
            StreamLine.SetActive(true);
        }
        else
        {
            StreamLine.SetActive(false);
        }
    }

    private void UpdateVelocity()
    {
        if (showVelocity)
        {
            Velocity.GetComponent<LineRenderer>().SetPosition(0, Vector3.zero);
            Velocity.GetComponent<LineRenderer>().SetPosition(1, data.velocity);
            Velocity.SetActive(true);
        }
        else
        {
            Velocity.SetActive(false);
        }
    }

    public void UpdateCell(CellData cellData)
    {
        data = cellData;
        data.newVelocity.Scale(gameObject.transform.localScale);
        data.velocity = data.newVelocity;
    }

    private void OnCollisionEnter(Collision collision)
    {
        data.permeability = 0;
    }
    private void OnCollisionExit(Collision collision)
    {
        data.permeability = 1;
    }

    private void OnCollisionStay(Collision collision)
    {
        data.permeability = 0;
    }

    private void OnTriggerEnter(Collider other)
    {
        data.permeability = 0;
    }

    private void OnTriggerExit(Collider other)
    {
           data.permeability = 1;
    }

    private void OnTriggerStay(Collider other)
    {
        data.permeability = 0;
    }
}
