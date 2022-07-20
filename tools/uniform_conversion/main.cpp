#include <iostream>
#include <vector>
#include <hdf5.h>
#include <mpi.h>
#include <string>
#include <sstream>
#include <cstdio>
#include <fstream>
#include <algorithm>    // std::sort
#include <fstream>
#include <filesystem>
#include <cmath>
namespace fs = std::filesystem;

#define BS 8
#define Cfactor 4

struct BlockGroup
{
  double h;
  int nx,ny,nz;
  double ox,oy,oz; //origin
  int level;
  int index[3];
};

void decompose_1D(int tasks, int & my_start, int & my_end, const int compression=1)
{
  int rank,size;
  MPI_Comm_rank(MPI_COMM_WORLD,&rank);
  MPI_Comm_size(MPI_COMM_WORLD,&size);
  if (tasks % compression != 0){
	  std::cout << "Tasks must be divisible by Cfactor.\n";
	  MPI_Abort(MPI_COMM_WORLD,1);
  }

  int my_share = tasks / compression / size;
  if (tasks % size != 0 && rank == size - 1) //last rank gets what's left
  {
    my_share += (tasks/compression) % size;
  }

  my_start = rank * (tasks/compression / size) * compression;
  my_end = my_start + my_share*compression;
}

std::vector<BlockGroup> get_amr_groups(std::string filename)
{
  hid_t file_id,  fapl_id;
  hid_t dataset_origins, fspace_origins;
  hid_t dataset_indices, fspace_indices;

  H5open();

  //open amr raw data - every rank will read that
  fapl_id = H5Pcreate(H5P_FILE_ACCESS);
  H5Pset_fapl_mpio(fapl_id, MPI_COMM_WORLD, MPI_INFO_NULL);
  file_id = H5Fopen((filename+"-groups.h5").c_str(), H5F_ACC_RDONLY, fapl_id);
  H5Pclose(fapl_id);

  dataset_origins = H5Dopen2(file_id, "origins", H5P_DEFAULT);
  dataset_indices = H5Dopen2(file_id, "indices", H5P_DEFAULT);

  //Read size of input dataset and store it to dim
  hsize_t dim_origins;
  hsize_t dim_indices;
  fspace_origins = H5Dget_space(dataset_origins);
  fspace_indices = H5Dget_space(dataset_indices);
  H5Sget_simple_extent_dims(fspace_origins, &dim_origins, NULL);
  H5Sget_simple_extent_dims(fspace_indices, &dim_indices, NULL);

  //Allocate vector to read dataset
  std::vector<double> origins(dim_origins);
  std::vector<int   > indices(dim_indices);
  H5Dread(dataset_origins, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, origins.data());
  H5Dread(dataset_indices, H5T_NATIVE_INT   , H5S_ALL, H5S_ALL, H5P_DEFAULT, indices.data());

  H5Dclose(dataset_origins);
  H5Dclose(dataset_indices);
  H5Sclose(fspace_origins);
  H5Sclose(fspace_indices);
  H5Fclose(file_id);
  H5close();

  std::vector<BlockGroup> groups(dim_origins/4);
  for (size_t i = 0 ; i < groups.size() ; i++)
  {
    groups[i].ox = origins[4*i  ]; 
    groups[i].oy = origins[4*i+1]; 
    groups[i].oz = origins[4*i+2]; 
    groups[i].h  = origins[4*i+3]; 
    groups[i].nx       = indices[7*i  ];
    groups[i].ny       = indices[7*i+1];
    groups[i].nz       = indices[7*i+2];
    groups[i].index[0] = indices[7*i+3];
    groups[i].index[1] = indices[7*i+4];
    groups[i].index[2] = indices[7*i+5];
    groups[i].level    = indices[7*i+6];
  }
  return groups;
}

std::vector<float> get_amr_dataset(std::string filename)
{
  hid_t file_id, dataset_id, fspace_id, fapl_id;

  H5open();

  //open amr raw data - every rank will read that
  fapl_id = H5Pcreate(H5P_FILE_ACCESS);
  H5Pset_fapl_mpio(fapl_id, MPI_COMM_WORLD, MPI_INFO_NULL);
  file_id = H5Fopen((filename+".h5").c_str(), H5F_ACC_RDONLY, fapl_id);
  H5Pclose(fapl_id);

  dataset_id = H5Dopen2(file_id, "dset", H5P_DEFAULT);

  //Read size of input dataset and store it to dim
  hsize_t dim;
  fspace_id = H5Dget_space(dataset_id);
  H5Sget_simple_extent_dims(fspace_id, &dim, NULL);

  //Allocate vector to read dataset
  std::vector<float> amr(dim);

  #if 1 //large datasets are read like this. HDF5 complains otherwise.
    const int max_chunk = 8;
    for (int chunk = 0 ; chunk < max_chunk ; chunk ++)
    {
      hsize_t count[1] = {dim/max_chunk};
      hsize_t base_tmp[1] = {chunk*count[0]};
  
      hid_t mspace_id = H5Screate_simple(1, count, NULL);
  
      H5Sselect_hyperslab(fspace_id, H5S_SELECT_SET, base_tmp, NULL, count, NULL);
  
      H5Dread(dataset_id, H5T_NATIVE_FLOAT, mspace_id, fspace_id, H5P_DEFAULT, amr.data()+base_tmp[0]);
  
      H5Sclose(mspace_id);
    }
  #else
    H5Dread(dataset_id, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT, amr.data());
  #endif

  H5Dclose(dataset_id);
  H5Sclose(fspace_id);
  H5Fclose(file_id);
  H5close();
  return amr;
}

void convert_to_uniform(std::string filename,int tttt)
{
  int rank,size;
  MPI_Comm_rank(MPI_COMM_WORLD,&rank);
  MPI_Comm_size(MPI_COMM_WORLD,&size);

  std::vector<float> amr = get_amr_dataset(filename);

  std::vector<BlockGroup> allGroups = get_amr_groups(filename);
  std::vector<int> base(allGroups.size());
  base[0] = 0;
  for (size_t i = 1 ; i < allGroups.size() ; i++)
  {
    base[i] = base[i-1] + allGroups[i-1].nx*allGroups[i-1].ny*allGroups[i-1].nz;
  }

  //find min h, domain extent and maxLevel
  double minh = 1e6;
  int    levelMax = -1;
  int    points[3] = {0,0,0};

  int my_start,my_end;
  decompose_1D(allGroups.size(),my_start,my_end);

  for (int i = my_start ; i < my_end ; i++)
  {
    minh = std::min(allGroups[i].h,minh);
    levelMax = std::max(allGroups[i].level,levelMax);
  }
  levelMax ++;
  MPI_Allreduce(MPI_IN_PLACE, &minh    , 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);
  MPI_Allreduce(MPI_IN_PLACE, &levelMax, 1, MPI_INT   , MPI_MAX, MPI_COMM_WORLD);

  for (int i = my_start ; i < my_end ; i++)
  {
    const int aux = 1 << (levelMax - 1 - allGroups[i].level);
    points[0] = std::max(points[0], (allGroups[i].index[0]*BS + allGroups[i].nx)*aux );
    points[1] = std::max(points[1], (allGroups[i].index[1]*BS + allGroups[i].ny)*aux );
    points[2] = std::max(points[2], (allGroups[i].index[2]*BS + allGroups[i].nz)*aux );
  }
  MPI_Allreduce(MPI_IN_PLACE, &points, 3, MPI_INT , MPI_MAX, MPI_COMM_WORLD);

  //the uniform domain is decomposed in the z-direction only!
  decompose_1D(points[2],my_start,my_end,Cfactor);

  size_t unc = (my_end-my_start)/Cfactor;
  unc *= points[1]*points[0]/Cfactor/Cfactor;
  std::vector<float> uniform_grid(unc);

  const double coef = 1.0 / Cfactor / Cfactor / Cfactor;

  #pragma omp parallel for
  for (size_t i = 0 ; i < allGroups.size() ; i++)
  {
    const BlockGroup & group = allGroups[i];
 
    const int aux = 1 << (levelMax - 1 - group.level);
    const int start_x = group.index[0]*BS*aux;
    const int start_y = group.index[1]*BS*aux;
    const int start_z = group.index[2]*BS*aux;

    //check if this group is within my part of the uniform domain
    const int end_z   = start_z + group.nz*aux;
    if (end_z < my_start || start_z > my_end) continue;

    for (int z = 0; z < group.nz; z++)
    for (int y = 0; y < group.ny; y++)
    for (int x = 0; x < group.nx; x++)
    {
      const float value = coef * amr[base[i] + x + y * group.nx + z*group.nx*group.ny];

      for (int z_up = aux * z; z_up < aux * (z+1); z_up++)
      for (int y_up = aux * y; y_up < aux * (y+1); y_up++)
      for (int x_up = aux * x; x_up < aux * (x+1); x_up++)
      {
        const int uniform_z = start_z + z_up;
        const int uniform_y = start_y + y_up;
        const int uniform_x = start_x + x_up;
        if (uniform_z >= my_start && uniform_z < my_end)
        {
          const int base = uniform_x/Cfactor + (uniform_y/Cfactor)*(points[0]/Cfactor) + ((uniform_z- my_start)/Cfactor)*(points[0]/Cfactor)*(points[1]/Cfactor);        
          uniform_grid[base] += value;
        }
      }
    }
  }

  if (rank == 0)
  {
    std::cout << "uniform domain size=" << points[0] /Cfactor << " x " << points[1] /Cfactor << " x " << points[2] /Cfactor << std::endl;

    std::stringstream s;
    s << "<?xml version=\"1.0\" ?>\n";
    s << "<!DOCTYPE Xdmf SYSTEM \"Xdmf.dtd\" []>\n";
    s << "<Xdmf Version=\"2.0\">\n";
    s << "<Domain>\n";
    s << "  <Time Value=\"" << std::scientific << 0.05*tttt << "\"/>\n\n";
    s << "  <Grid GridType=\"Uniform\">\n";
    s << "    <Topology TopologyType=\"3DCoRectMesh\" Dimensions=\" " << points[2]/Cfactor + 1 << " " << points[1]/Cfactor + 1<< " " << points[0]/Cfactor + 1 << "\"/>\n";
    s << "    <Geometry GeometryType=\"ORIGIN_DXDYDZ\">\n";
    s << "       <DataItem Dimensions=\"3\" NumberType=\"Double\" Precision=\"8\" " "Format=\"XML\">\n";
    s << "            " << std::scientific << 0.0 << " " << 0.0 << " " << 0.0 << "\n";
    s << "       </DataItem>\n";
    s << "      <DataItem Dimensions=\"3\" NumberType=\"Double\" Precision=\"8\" " "Format=\"XML\">\n";
    s << "            " << std::scientific << Cfactor*minh <<" "<< Cfactor*minh <<" "<< Cfactor*minh << "\n";
    s << "       </DataItem>\n";
    s << "   </Geometry>\n";
    s << "   <Attribute Name=\"data\" AttributeType=\"" << "Scalar"<< "\" Center=\"Cell\">\n";
    s << "      <DataItem ItemType=\"Uniform\"  Dimensions=\" " << points[2]/Cfactor << " " << points[1]/Cfactor << " " << points[0]/Cfactor << " " << "\" NumberType=\"Float\" Precision=\" " << (int)sizeof(H5T_NATIVE_FLOAT) << "\" Format=\"HDF\">\n";
    s << "       " << (filename + "-uniform.h5").c_str() << ":/" << "data" << "\n";
    s << "     </DataItem>\n";
    s << "   </Attribute>\n";  
    s << "  </Grid>\n\n";
    s << "</Domain>\n";
    s << "</Xdmf>\n";
    std::string st = s.str();
    std::ofstream out((filename + "-uniform.xmf").c_str());
    out << st;
    out.close();
  }

  //dump uniform grid
  {
    hid_t file_id, dataset_id, fspace_id, fapl_id, mspace_id;
    H5open();
    fapl_id = H5Pcreate(H5P_FILE_ACCESS);
    H5Pset_fapl_mpio(fapl_id, MPI_COMM_WORLD, MPI_INFO_NULL);
    file_id = H5Fcreate((filename+"-uniform.h5").c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, fapl_id);        
    H5Pclose(fapl_id);
    fapl_id = H5Pcreate(H5P_DATASET_XFER);

    H5Pset_dxpl_mpio(fapl_id, H5FD_MPIO_COLLECTIVE);
    hsize_t dims[3]  = { (hsize_t)points[2]/Cfactor,(hsize_t)points[1]/Cfactor,(hsize_t)points[0]/Cfactor };

    ////compressed dataset
    hid_t plist_id = H5Pcreate(H5P_DATASET_CREATE);
    #if 0
        hsize_t cdims[3];
        cdims[0] = dims[0] / 8;//64;
        cdims[1] = dims[1] / 8;//64;
        cdims[2] = dims[2] / 8;//64;
        if (dims[0] % 64 == 0 && dims[1] % 64 == 0 && dims[2] % 64 == 0)
        {
          H5Pset_chunk(plist_id, 3, cdims);
          H5Pset_deflate(plist_id, 5);
          if (rank == 0)
            std::cout << " -> data compression enabled." << std::endl;
        }
        else
        {
          if (rank == 0)
            std::cout << " -> data compression disabled." << std::endl;
        }
    #endif

    fspace_id        = H5Screate_simple(3, dims, NULL);
    dataset_id       = H5Dcreate (file_id, "data", H5T_NATIVE_FLOAT ,fspace_id,H5P_DEFAULT,plist_id,H5P_DEFAULT);
    H5Sclose(fspace_id);

    fspace_id = H5Dget_space(dataset_id);

    hsize_t count[3] = {(hsize_t)(my_end-my_start)/Cfactor,(hsize_t)points[1]/Cfactor,(hsize_t)points[0]/Cfactor};
    hsize_t base_tmp[3] = {(hsize_t)my_start/Cfactor,0,0};

    mspace_id = H5Screate_simple(3, count, NULL);
    H5Sselect_hyperslab(fspace_id, H5S_SELECT_SET, base_tmp, NULL, count, NULL);
    H5Dwrite(dataset_id, H5T_NATIVE_FLOAT,mspace_id,fspace_id,fapl_id,uniform_grid.data());

    H5Sclose(mspace_id);
    H5Sclose(fspace_id);
    H5Dclose(dataset_id);
    H5Pclose(fapl_id);
    H5Fclose(file_id);
    H5Pclose(plist_id);
    H5close();
  } 
}


int main(int argc, char **argv)
{
  int provided;
  const auto SECURITY = MPI_THREAD_FUNNELED;
  MPI_Init_thread(&argc, &argv, SECURITY, &provided);
  if (provided < SECURITY ) {
    printf("ERROR: MPI implementation does not have required thread support\n");
    fflush(0); MPI_Abort(MPI_COMM_WORLD, 1);
  }
  int rank;
  MPI_Comm_rank(MPI_COMM_WORLD,&rank);

  std::vector<std::string> filenames;
  std::string path("./");
  std::string ext(".h5");
  for (auto &p : fs::recursive_directory_iterator(path))
  {
    if (p.path().extension() == ext)
    {
        std::string s = p.path().stem().string();
        if ( s.back() != 's' && s.back() != 'm')
        {
          filenames.push_back(p.path().stem().string());
        }
    }
  }
  std::sort(filenames.begin(),filenames.end());

  for (int i = filenames.size()-1 ; i >= 0  ; i --)
  {
    if (rank == 0)
	    std::cout << "processing files: " << filenames[i] << std::endl;
    convert_to_uniform(filenames[i],i);
    MPI_Barrier(MPI_COMM_WORLD);
  }
  MPI_Finalize();
  return 0;
}