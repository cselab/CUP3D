#include <iostream>
#include <vector>
#include <hdf5.h>
#include <mpi.h>
#include <string>
#include <sstream>
#include <cstdio>
#include <fstream>

struct BlockGroup
{
  double h;
  int nx,ny,nz;
  double ox,oy,oz; //origin
  int level;
  int index[3];
};

void decompose_1D(long long tasks, long long & my_start, long long & my_end)
{
  int rank,size;
  MPI_Comm_rank(MPI_COMM_WORLD,&rank);
  MPI_Comm_size(MPI_COMM_WORLD,&size);

  long long my_share = tasks / size;
  if (tasks % size != 0 && rank == size - 1) //last rank gets what's left
  {
    my_share += tasks % size;
  }

  my_start = rank * (tasks / size);
  my_end = my_start + my_share;
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

std::vector<double> get_amr_dataset(std::string filename)
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
  std::vector<double> amr(dim);
  H5Dread(dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, amr.data());

  H5Dclose(dataset_id);
  H5Sclose(fspace_id);
  H5Fclose(file_id);
  H5close();
  return amr;
}

void convert_to_uniform(std::string filename)
{
  const int BS = 8;

  int rank,size;
  MPI_Comm_rank(MPI_COMM_WORLD,&rank);
  MPI_Comm_size(MPI_COMM_WORLD,&size);

  std::vector<double> amr = get_amr_dataset(filename);

  std::vector<BlockGroup> allGroups = get_amr_groups(filename);
  std::vector<long long> base(allGroups.size());
  base[0] = 0;
  for (size_t i = 1 ; i < allGroups.size() ; i++)
  {
    base[i] = base[i-1] + allGroups[i-1].nx*allGroups[i-1].ny*allGroups[i-1].nz;
  }

  //find min h, domain extent and maxLevel
  double minh = 1e6;
  int    levelMax = -1;
  long long points[3] = {0,0,0};

  long long my_start,my_end;
  decompose_1D(allGroups.size(),my_start,my_end);

  for (long long i = my_start ; i < my_end ; i++)
  {
    minh = std::min(allGroups[i].h,minh);
    levelMax = std::max(allGroups[i].level,levelMax);
  }
  levelMax ++;
  MPI_Allreduce(MPI_IN_PLACE, &minh    , 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);
  MPI_Allreduce(MPI_IN_PLACE, &levelMax, 1, MPI_INT   , MPI_MAX, MPI_COMM_WORLD);

  for (long long i = my_start ; i < my_end ; i++)
  {
    const int aux = 1 << (levelMax - 1 - allGroups[i].level);
    points[0] = std::max(points[0], (long long)(allGroups[i].index[0]*BS + allGroups[i].nx)*aux );
    points[1] = std::max(points[1], (long long)(allGroups[i].index[1]*BS + allGroups[i].ny)*aux );
    points[2] = std::max(points[2], (long long)(allGroups[i].index[2]*BS + allGroups[i].nz)*aux );
  }
  MPI_Allreduce(MPI_IN_PLACE, &points, 3, MPI_LONG_LONG , MPI_MAX, MPI_COMM_WORLD);

  if (rank == 0)
    std::cout << "uniform domain size=" << points[0] << " x " << points[1] << " x " << points[2] << std::endl;

  //the uniform domain is decomposed in the x-direction only!
  decompose_1D(points[0],my_start,my_end);
  std::vector<float> uniform_grid((my_end-my_start)*points[1]*points[2]);

  #pragma omp parallel for
  for (size_t i = 0 ; i < allGroups.size() ; i++)
  {
    const BlockGroup & group = allGroups[i];
    //check if this group is within my part of the uniform domain

    const int aux = 1 << (levelMax - 1 - group.level);
    const long long start_x = group.index[0]*BS*aux;
    const long long end_x   = start_x + group.nx*aux;

    if (end_x < my_start) continue;
    if (start_x > my_end) continue;

    const long long start_y = group.index[1]*BS*aux;
    const long long start_z = group.index[2]*BS*aux;


    for (int z = 0; z < group.nz; z++)
    for (int y = 0; y < group.ny; y++)
    for (int x = 0; x < group.nx; x++)
    {
      const double value = amr[base[i] + x + y * group.nx + z*group.nx*group.ny];

      for (int z_up = aux * z; z_up < aux * (z+1); z_up++)
      for (int y_up = aux * y; y_up < aux * (y+1); y_up++)
      for (int x_up = aux * x; x_up < aux * (x+1); x_up++)
      {
        const long long uniform_x = start_x + x_up;
        if (uniform_x >= my_start && uniform_x < my_end)
        {
          const long long uniform_y = start_y + y_up;
          const long long uniform_z = start_z + z_up;
          const int base_up = uniform_x - my_start + uniform_y*(my_end-my_start) + uniform_z*(my_end-my_start)*points[1];        
          uniform_grid[base_up] = (float)value;
        }
      }
    }
  }
  if (rank == 0)
  {
      std::stringstream s;
      s << "<?xml version=\"1.0\" ?>\n";
      s << "<!DOCTYPE Xdmf SYSTEM \"Xdmf.dtd\" []>\n";
      s << "<Xdmf Version=\"2.0\">\n";
      s << "<Domain>\n";
      //s << "  <Time Value=\"" << std::scientific << absTime << "\"/>\n\n";
      s << "  <Grid GridType=\"Uniform\">\n";
      s << "    <Topology TopologyType=\"3DCoRectMesh\" Dimensions=\" " << points[2] + 1 << " " << points[1] + 1<< " " << points[0] + 1 << "\"/>\n";
      s << "    <Geometry GeometryType=\"ORIGIN_DXDYDZ\">\n";
      s << "       <DataItem Dimensions=\"3\" NumberType=\"Double\" Precision=\"8\" " "Format=\"XML\">\n";
      s << "            " << std::scientific << 0.0 << " " << 0.0 << " " << 0.0 << "\n";
      s << "       </DataItem>\n";
      s << "      <DataItem Dimensions=\"3\" NumberType=\"Double\" Precision=\"8\" " "Format=\"XML\">\n";
      s << "            " << std::scientific << minh <<" "<< minh <<" "<< minh << "\n";
      s << "       </DataItem>\n";
      s << "   </Geometry>\n";
      s << "   <Attribute Name=\"data\" AttributeType=\"" << "Scalar"<< "\" Center=\"Cell\">\n";
      s << "      <DataItem ItemType=\"Uniform\"  Dimensions=\" " << points[2] << " " << points[1] << " " << points[0] << " " << "\" NumberType=\"Float\" Precision=\" " << (int)sizeof(H5T_NATIVE_FLOAT) << "\" Format=\"HDF\">\n";
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

    //hsize_t chunk_dims[3];
    //chunk_dims[0] = 8;   
    //chunk_dims[1] = 8;   
    //chunk_dims[2] = 8;   
    //hid_t memspace  = H5Screate_simple(3, chunk_dims, NULL); 
    //H5Pset_chunk(fapl_id, 3, chunk_dims);
    //H5Sclose(memspace);

    H5Pset_dxpl_mpio(fapl_id, H5FD_MPIO_COLLECTIVE);
    hsize_t dims[3]  = { (hsize_t)points[2],(hsize_t)points[1],(hsize_t)points[0] };
    fspace_id        = H5Screate_simple(3, dims, NULL);
    dataset_id       = H5Dcreate (file_id, "data", H5T_NATIVE_FLOAT ,fspace_id,H5P_DEFAULT,fapl_id,H5P_DEFAULT);
    H5Sclose(fspace_id);

    fspace_id = H5Dget_space(dataset_id);

    hsize_t count[3] = {(hsize_t)points[2],(hsize_t)points[1],(hsize_t)(my_end-my_start)};
    hsize_t base_tmp[3] = {0,0,(hsize_t)my_start};
    mspace_id = H5Screate_simple(3, count, NULL);
    H5Sselect_hyperslab(fspace_id, H5S_SELECT_SET, base_tmp, NULL, count, NULL);
    H5Dwrite(dataset_id, H5T_NATIVE_FLOAT,mspace_id,fspace_id,fapl_id,uniform_grid.data());
    H5Sclose(mspace_id);

    H5Sclose(fspace_id);
    H5Dclose(dataset_id);
    H5Pclose(fapl_id);
    H5Fclose(file_id);
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
  filenames.push_back("chi_init000000000");

  for (size_t i = 0 ; i < filenames.size() ; i ++)
  {
    std::cout << "processing file: " << filenames[i] << std::endl;
    convert_to_uniform(filenames[i]);
    MPI_Barrier(MPI_COMM_WORLD);
  }
  MPI_Finalize();
  return 0;
}