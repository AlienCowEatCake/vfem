#include "vfem.h"

// Вывод данных в 3D сетке
bool VFEM::output(const string & tecplot_filename)
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
        cout << "Error in " << __FILE__ << ":" << __LINE__
             << " while writing file " << tecplot_filename << endl;
        return false;
    }

    tecplot_file.precision(17);
    tecplot_file.setf(ios::scientific);

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
            tecplot_file << fes[i].nodes[j]->num + 1 << " ";
        tecplot_file << "\n";
    }

    tecplot_file << "\n";
    tecplot_file.flush();
    tecplot_file.close();
    return true;
}

// Вывод данных в 2D сетке
bool VFEM::output_slice(const string & tecplot_filename, char slice_var, double slice_val,
                        char var1, double min_var1, double max_var1, size_t num_var_1,
                        char var2, double min_var2, double max_var2, size_t num_var_2)
{
    if(min_var1 > max_var1) swap(max_var1, min_var1);
    if(min_var2 > max_var2) swap(max_var2, min_var2);

    double step_var1 = (max_var1 - min_var1) / (double)(num_var_1 - 1);
    double step_var2 = (max_var2 - min_var2) / (double)(num_var_2 - 1);
    size_t index_slice = 0, index_1 = 0, index_2 = 0;

    // Определяем по какой переменной сечение
    if      (slice_var == 'x' || slice_var == 'X') index_slice = 0;
    else if (slice_var == 'y' || slice_var == 'Y') index_slice = 1;
    else if (slice_var == 'z' || slice_var == 'Z') index_slice = 2;
    else
    {
        cout << "Unknown variable, breaking ..." << endl;
        return false;
    }

    // Определяем, какая переменная первая
    if      (var1 == 'x' || var1 == 'X') index_1 = 0;
    else if (var1 == 'y' || var1 == 'Y') index_1 = 1;
    else if (var1 == 'z' || var1 == 'Z') index_1 = 2;
    else
    {
        cout << "Unknown variable, breaking ..." << endl;
        return false;
    }

    // Определяем, какая переменная вторая
    if      (var2 == 'x' || var2 == 'X') index_2 = 0;
    else if (var2 == 'y' || var2 == 'Y') index_2 = 1;
    else if (var2 == 'z' || var2 == 'Z') index_2 = 2;
    else
    {
        cout << "Unknown variable, breaking ..." << endl;
        return false;
    }

    cout << "Writing slice to Tecplot ..." << endl;

    ofstream tecplot_file;
    tecplot_file.open(tecplot_filename.c_str(), ios::out);

    if(!tecplot_file.good())
    {
        cout << "Error in " << __FILE__ << ":" << __LINE__
             << " while writing file " << tecplot_filename << endl;
        return false;
    }

    tecplot_file << "TITLE = \"Slice " << slice_var << " = " << slice_val << "\"\n";
    tecplot_file << "VARIABLES = \"" << var1 <<"\", \"" << var2 << "\", \"ExR\", \"EyR\", \"EzR\", \"ExI\", \"EyI\", \"EzI\", \"abs(E)\", \"div(ER)\", \"div(EI)\"\n";
    tecplot_file << "ZONE I= " << num_var_2 << ", J= " << num_var_1 << ", F=POINT\n";

    tecplot_file.precision(17);
    tecplot_file.setf(ios::scientific);

    point p(0, 0, 0);
    p[index_slice] = slice_val;
    for(size_t i = 0; i < num_var_1; i++)
    {
        double v1 = min_var1 + step_var1 * (double)i;
        for(size_t j = 0; j < num_var_2; j++)
        {
            double v2 = min_var2 + step_var2 * (double)j;
            p[index_1] = v1;
            p[index_2] = v2;
            const finite_element * fe = get_fe(p);
            cvector3 sol = solution(p, fe);
            complex<double> div_sol = divergence(p, fe) * complex<double>(0, 1) / fe->phys->omega;
            tecplot_file << v1 << " " << v2 << " "
                         << sol.x.real() << " " << sol.y.real() << " " << sol.z.real() << " "
                         << sol.x.imag() << " " << sol.y.imag() << " " << sol.z.imag() << " "
                         << sol.norm() << " " << div_sol.real() << " " << div_sol.imag() << "\n";
        }
    }

    tecplot_file << "\n";
    tecplot_file.flush();
    tecplot_file.close();
    return true;
}

// Вывод данных по линии
bool VFEM::output_line(const string & tecplot_filename, char line_var1, double line_val1,
                       char line_var2, double line_val2, char var3, double min_var3,
                       double max_var3, size_t num_var)
{
    if(min_var3 > max_var3) swap(max_var3, min_var3);
    double step_var3 = (max_var3 - min_var3) / (double)(num_var - 1);
    size_t index1 = 0, index2 = 0, index_3 = 0;

    if      (line_var1 == 'x' || line_var1 == 'X') index1 = 0;
    else if (line_var1 == 'y' || line_var1 == 'Y') index1 = 1;
    else if (line_var1 == 'z' || line_var1 == 'Z') index1 = 2;
    else
    {
        cout << "Unknown variable, breaking ..." << endl;
        return false;
    }

    if      (line_var2 == 'x' || line_var2 == 'X') index2 = 0;
    else if (line_var2 == 'y' || line_var2 == 'Y') index2 = 1;
    else if (line_var2 == 'z' || line_var2 == 'Z') index2 = 2;
    else
    {
        cout << "Unknown variable, breaking ..." << endl;
        return false;
    }

    if      (var3 == 'x' || var3 == 'X') index_3 = 0;
    else if (var3 == 'y' || var3 == 'Y') index_3 = 1;
    else if (var3 == 'z' || var3 == 'Z') index_3 = 2;
    else
    {
        cout << "Unknown variable, breaking ..." << endl;
        return false;
    }

    cout << "Writing line to Tecplot ..." << endl;

    ofstream tecplot_file;
    tecplot_file.open(tecplot_filename.c_str(), ios::out);

    if(!tecplot_file.good())
    {
        cout << "Error in " << __FILE__ << ":" << __LINE__
             << " while writing file " << tecplot_filename << endl;
        return false;
    }

    tecplot_file << "TITLE = \"Line " << line_var1 << " = " << line_val1 << ", " << line_var2 << " = " << line_val2 << "\"\n";
    tecplot_file << "VARIABLES = \"" << var3 << "\", \"ExR\", \"EyR\", \"EzR\", \"ExI\", \"EyI\", \"EzI\", \"abs(E)\"\n";
    tecplot_file << "ZONE I= " << num_var << ", F=POINT\n";

    tecplot_file.precision(17);
    tecplot_file.setf(ios::scientific);

    point p(0, 0, 0);
    p[index1] = line_val1;
    p[index2] = line_val2;
    for(size_t i = 0; i < num_var; i++)
    {
        double v3 = min_var3 + step_var3 * (double)i;
        p[index_3] = v3;
        cvector3 sol = solution(p);
        tecplot_file << v3 << " "
                   << sol.x.real() << " " << sol.y.real() << " " << sol.z.real() << " "
                   << sol.x.imag() << " " << sol.y.imag() << " " << sol.z.imag() << " "
                   << sol.norm() << "\n";
    }

    tecplot_file << "\n";
    tecplot_file.flush();
    tecplot_file.close();
    return true;
}
