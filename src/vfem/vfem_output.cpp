#include "vfem.h"

void VFEM::output(const string & tecplot_filename)
{
    cout << "Writing to Tecplot ..." << endl;

    size_t fes_num = fes.size();
    vector<cvector3> E_vals;
    E_vals.resize(fes_num);
    for(size_t k = 0; k < fes_num; k++)
        E_vals[k] = solution(fes[k].barycenter, &fes[k]);

    ofstream tecplot_file;
    tecplot_file.open(tecplot_filename.c_str(), ios::out);

    if(!tecplot_file.good())
    {
        cerr << "Error in " << __FILE__ << ":" << __LINE__
             << " while writing file " << tecplot_filename << endl;
        throw IO_FILE_ERROR;
    }

    tecplot_file.precision(17);
    tecplot_file.setf(ios::fixed);

    size_t nodes_num = nodes.size();
    tecplot_file << "TITLE     = \"" << "Title" << "\" \n  VARIABLES = \"x\" \n \"y\" \n \"z\" \n \"ExR\" \n \"EyR\" \n \"EzR\" \n \"ExI\" \n \"EyI\" \n \"EzI\" \n \"abs(E)\"\n ";
    tecplot_file << "ZONE T=\"ZONE1\" \n N = " << nodes_num << " E = " << fes_num << " , F=FEBLOCK ET=Tetrahedron \nVARLOCATION = (NODAL NODAL NODAL CELLCENTERED CELLCENTERED CELLCENTERED CELLCENTERED CELLCENTERED CELLCENTERED CELLCENTERED)\n";

    for(size_t i = 0; i < nodes_num; i++)
        tecplot_file << nodes[i].x << "\n";
    for(size_t i = 0; i < nodes_num; i++)
        tecplot_file << nodes[i].y << "\n";
    for(size_t i = 0; i < nodes_num; i++)
        tecplot_file << nodes[i].z << "\n";

    for(size_t i = 0; i < fes_num; i++)
        tecplot_file << E_vals[i].x.real() << "\n";
    for(size_t i = 0; i < fes_num; i++)
        tecplot_file << E_vals[i].y.real() << "\n";
    for(size_t i = 0; i < fes_num; i++)
        tecplot_file << E_vals[i].z.real() << "\n";

    for(size_t i = 0; i < fes_num; i++)
        tecplot_file << E_vals[i].x.imag() << "\n";
    for(size_t i = 0; i < fes_num; i++)
        tecplot_file << E_vals[i].y.imag() << "\n";
    for(size_t i = 0; i < fes_num; i++)
        tecplot_file << E_vals[i].z.imag() << "\n";

    for(size_t i = 0; i < fes_num; i++)
        tecplot_file << E_vals[i].norm() << "\n";

    for(size_t i = 0; i < fes_num; i++)
    {
        for(size_t j = 0; j < 4; j++)
            tecplot_file << fes[i].get_node(j).num + 1 << " ";
        tecplot_file << "\n";
    }

    tecplot_file << "\n";
    tecplot_file.flush();
    tecplot_file.close();
}

void VFEM::output_slice(const string & tecplot_filename, char slice_var, double slice_val,
                        char var1, double min_var1, double max_var1, double step_var1,
                        char var2, double min_var2, double max_var2, double step_var2)
{
    if(min_var1 > max_var1) swap(max_var1, min_var1);
    if(min_var2 > max_var2) swap(max_var2, min_var2);

    double min_x = 0.0, max_x = 0.0, step_x = 0.0;
    double min_y = 0.0, max_y = 0.0, step_y = 0.0;
    double min_z = 0.0, max_z = 0.0, step_z = 0.0;
    size_t nx = 0, ny = 0, nz = 0;

    if(slice_var == 'x' || slice_var == 'X')
    {
        min_x = max_x = slice_val;
        nx = 0;
        step_x = 0.0;
    }
    else if (slice_var == 'y' || slice_var == 'Y')
    {
        min_y = max_y = slice_val;
        ny = 0;
        step_y = 0.0;
    }
    else if (slice_var == 'z' || slice_var == 'Z')
    {
        min_z = max_z = slice_val;
        nz = 0;
        step_z = 0.0;
    }
    else
    {
        cerr << "Unknown variable, breaking ..." << endl;
        return;
    }

    if(var1 == 'x' || var1 == 'X')
    {
        min_x = min_var1;
        max_x = max_var1;
        step_x = step_var1;
        nx = (size_t)((max_var1 - min_var1) / step_var1);
    }
    else if (var1 == 'y' || var1 == 'Y')
    {
        min_y = min_var1;
        max_y = max_var1;
        step_y = step_var1;
        ny = (size_t)((max_var1 - min_var1) / step_var1);
    }
    else if (var1 == 'z' || var1 == 'Z')
    {
        min_z = min_var1;
        max_z = max_var1;
        step_z = step_var1;
        nz = (size_t)((max_var1 - min_var1) / step_var1);
    }
    else
    {
        cerr << "Unknown variable, breaking ..." << endl;
        return;
    }

    if(var2 == 'x' || var2 == 'X')
    {
        min_x = min_var2;
        max_x = max_var2;
        step_x = step_var2;
        nx = (size_t)((max_var2 - min_var2) / step_var2);
    }
    else if (var2 == 'y' || var2 == 'Y')
    {
        min_y = min_var2;
        max_y = max_var2;
        step_y = step_var2;
        ny = (size_t)((max_var2 - min_var2) / step_var2);
    }
    else if (var2 == 'z' || var2 == 'Z')
    {
        min_z = min_var2;
        max_z = max_var2;
        step_z = step_var2;
        nz = (size_t)((max_var2 - min_var2) / step_var2);
    }
    else
    {
        cerr << "Unknown variable, breaking ..." << endl;
        return;
    }

    cout << "Writing slice to Tecplot ..." << endl;

    ofstream tecplot_file;
    tecplot_file.open(tecplot_filename.c_str(), ios::out);

    if(!tecplot_file.good())
    {
        cerr << "Error in " << __FILE__ << ":" << __LINE__
             << " while writing file " << tecplot_filename << endl;
        throw IO_FILE_ERROR;
    }

    tecplot_file << "TITLE = \"Slice " << slice_var << " = " << slice_val << "\"\n";
    tecplot_file << "VARIABLES = \"x\", \"y\", \"z\", \"ExR\", \"EyR\", \"EzR\", \"ExI\", \"EyI\", \"EzI\", \"abs(E)\"\n";
    tecplot_file << "ZONE I= " << nx + 1 << ", J= " << ny + 1 << ", K= " << nz + 1 << ", F=POINT\n";

    tecplot_file.precision(17);
    tecplot_file.setf(ios::scientific);

    for(size_t i = 0; i <= nz; i++)
    {
        double z = min_z + step_z * (double)i;
        for(size_t j = 0; j <= ny; j++)
        {
            double y = min_y + step_y * (double)j;
            for(size_t k = 0; k <= nx; k++)
            {
                double x = min_x + step_x * (double)k;
                cvector3 sol = solution(point(x, y, z));
                tecplot_file << x << " " << y << " " << z << " "
                             << sol.x.real() << " " << sol.y.real() << " " << sol.z.real() << " "
                             << sol.x.imag() << " " << sol.y.imag() << " " << sol.z.imag() << " "
                             << sol.norm() << "\n";
            }
        }
    }

    tecplot_file << "\n";
    tecplot_file.flush();
    tecplot_file.close();
}

