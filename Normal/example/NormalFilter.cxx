#include <itkVector.h>
#include <itkQuadEdgeMesh.h>
#include <itkVTKPolyDataReader.h>

#include "itkQuadEdgeMeshExtendedTraits.h"
#include "itkQENormalFilter.h"

using namespace itk;
using namespace std;

int main( int argc, char** argv )
{
    if( argc < 3 )
    {
      std::cout <<"***     Normal Filter     ***" <<std::endl;
      std::cout <<"It requires 2 arguments:" <<std::endl;
      std::cout <<" 1- Input mesh file name" <<std::endl;
      std::cout <<" 2- Weight Type" <<std::endl;
      std::cout <<"    * 0: GOURAUD" <<std::endl;
      std::cout <<"    * 1: THRUMER" <<std::endl;
      std::cout <<"    * 2: AREA" <<std::endl;
      std::cout <<"*****************************" <<std::endl;
      return EXIT_FAILURE;
    }
    const unsigned int Dimension = 3;
    typedef double CoordType;
    typedef QuadEdgeMesh< CoordType, Dimension > InputMeshType;

    typedef Vector< CoordType, Dimension > VectorType;

    typedef QuadEdgeMeshExtendedTraits <
            VectorType,
            Dimension,
            2,
            CoordType,
            CoordType,
            VectorType,
            bool,
            bool > Traits;
    typedef QuadEdgeMesh < VectorType, Dimension, Traits > OutputMeshType;
    //typedef OutputMeshType InputMeshType;

    typedef VTKPolyDataReader< InputMeshType > ReaderType;

    ReaderType::Pointer reader = ReaderType::New( );
    reader->SetFileName( argv[1] );
    try
    {
      reader->Update( );
    }
    catch( itk::ExceptionObject & exp )
    {
      std::cerr << "Exception thrown while reading the input file " << std::endl;
      std::cerr << exp << std::endl;
      return EXIT_FAILURE;
    }

    InputMeshType::Pointer mesh = reader->GetOutput( );

    typedef QENormalFilter< InputMeshType, OutputMeshType > NormalFilterType;
//    typedef QENormalFilter< OutputMeshType, OutputMeshType > NormalFilterType;
    NormalFilterType::Pointer normals = NormalFilterType::New( );
    normals->SetInput( mesh );

    NormalFilterType::WeightType weight_type;

    switch( atoi( argv[2] ) )
    {
      case 0:
        weight_type = NormalFilterType::GOURAUD;
        break;
      case 1:
        weight_type = NormalFilterType::THURMER;
        break;
      case 2:
        weight_type = NormalFilterType::AREA;
        break;
      default:
        return EXIT_FAILURE;
    }
    normals->SetWeight( weight_type );
    normals->Update( );

    OutputMeshType::Pointer output = normals->GetOutput( );

    OutputMeshType::PointDataContainerPointer pointdata = output->GetPointData( );

    cout <<"*********************************" <<endl;
    cout <<"Vertex Normal" <<endl;
    for( OutputMeshType::PointDataContainerIterator d_it = pointdata->Begin( );
         d_it != pointdata->End( );
         d_it++ )
    {
        cout <<d_it->Index( ) <<"  " <<d_it->Value( ) <<endl;
    }

    cout <<endl;
    cout <<"*********************************" <<endl;
    cout <<"Face Normal" <<endl;

    OutputMeshType::CellDataContainerPointer celldata = output->GetCellData( );


    for( OutputMeshType::CellDataContainerIterator n_it = celldata->Begin( );
         n_it != celldata->End( );
         n_it++ )
    {
        cout <<n_it->Index( ) <<"  " <<n_it->Value( ) <<endl;
    }

    return EXIT_SUCCESS;
}
