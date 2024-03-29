#ifndef __itkQENormalFilter_txx
#define __itkQENormalFilter_txx

#include "itkQENormalFilter.h"

namespace itk
{
  template< class TInputMesh, class TOutputMesh >
  QENormalFilter< TInputMesh, TOutputMesh >::
  QENormalFilter( ) : Superclass( ), m_Weight( THURMER )
  {}

  template< class TInputMesh, class TOutputMesh >
  QENormalFilter< TInputMesh, TOutputMesh >::
  ~QENormalFilter( )
  {}

  template< class TInputMesh, class TOutputMesh >
  typename QENormalFilter< TInputMesh, TOutputMesh >::OutputFaceNormalType
  QENormalFilter< TInputMesh, TOutputMesh >::
  ComputeFaceNormal( OutputPolygonType* iPoly )
  {
    OutputMeshPointer output = this->GetOutput( );

    OutputPointType pt[3];
    int k( 0 );

    OutputQEType* edge = iPoly->GetEdgeRingEntry( );
    OutputQEType* temp = edge;
    do
    {
      pt[k++] = output->GetPoint( temp->GetOrigin( ) );
      temp = temp->GetLnext( );
    }while( temp != edge );

    return TriangleType::ComputeNormal( pt[0], pt[1], pt[2] );
  }

  template< class TInputMesh, class TOutputMesh >
  void
  QENormalFilter< TInputMesh, TOutputMesh >::
  ComputeAllFaceNormals( )
  {
    OutputMeshPointer output = this->GetOutput( );
    OutputPolygonType* poly;

    for( OutputCellsContainerConstIterator
          cell_it = output->GetCells( )->Begin( );
         cell_it != output->GetCells( )->End( );
         ++cell_it )
    {
      poly = dynamic_cast< OutputPolygonType* >( cell_it.Value( ) );

      if( poly != 0 )
        if( poly->GetNumberOfPoints( ) == 3 )
          output->SetCellData( cell_it->Index( ), ComputeFaceNormal( poly ) );
    }
  }

  template< class TInputMesh, class TOutputMesh >
  void
  QENormalFilter< TInputMesh, TOutputMesh >::
  ComputeAllVertexNormals( )
  {
    OutputMeshPointer output = this->GetOutput( );
    OutputPointsContainerPointer points = output->GetPoints( );
    OutputPointIdentifier id;

    for( OutputPointsContainerIterator it = points->Begin( );
          it != points->End( );
         ++it )
    {
      id = it->Index( );
      output->SetPointData( id, ComputeVertexNormal( id ) );
    }
  }

  template< class TInputMesh, class TOutputMesh >
  typename QENormalFilter< TInputMesh, TOutputMesh >::OutputVertexNormalType
  QENormalFilter< TInputMesh, TOutputMesh >::
  ComputeVertexNormal( const OutputPointIdentifier& iId )
  {
    OutputMeshPointer output = this->GetOutput( );

    OutputQEType* edge = output->FindEdge( iId );
    OutputQEType* temp = edge;
    OutputCellIdentifier cell_id( 0 );

    OutputVertexNormalType n( 0. );
    OutputFaceNormalType face_normal( 0. );

    do
    {
      cell_id = temp->GetLeft( );
      output->GetCellData( cell_id, &face_normal );
      n += face_normal * Weight( iId, cell_id );
      temp = temp->GetOnext( );
    } while( temp != edge );

    n.Normalize( );
    return n;
  }

  template< class TInputMesh, class TOutputMesh >
  typename QENormalFilter< TInputMesh, TOutputMesh >::
    OutputVertexNormalComponentType
  QENormalFilter< TInputMesh, TOutputMesh >::
  Weight( const OutputPointIdentifier& iPId,
          const OutputCellIdentifier& iCId )
  {
    if( m_Weight == GOURAUD )
      return static_cast< OutputVertexNormalComponentType >( 1. );
    else
    {
      OutputMeshPointer output = this->GetOutput( );

      OutputPolygonType* poly = dynamic_cast< OutputPolygonType* >(
        output->GetCells( )->GetElement( iCId ) );
      if( poly != 0 ) // this test should be removed...
      {
        // this test should be removed...
        if( poly->GetNumberOfPoints( ) == 3 )
        {
          int internal_id( 0 ), k( 0 );
          OutputPointType pt[3];
          OutputVectorType u, v;

          OutputQEType* edge = poly->GetEdgeRingEntry( );
          OutputQEType* temp = edge;
          k = 0;
          do
          {
            pt[k] = output->GetPoint( temp->GetOrigin( ) );
            if( temp->GetOrigin( ) == iPId )
              internal_id = k;

            temp = temp->GetLnext( );
            k++;
          }while( temp != edge );

          switch( m_Weight )
          {
            default:
            case GOURAUD:
              return static_cast< OutputVertexNormalComponentType >( 1. );
            case THURMER:
              // this implementation may be included inside itkTriangle
              switch( internal_id )
              {
                case 0:
                  u = pt[1] - pt[0];
                  v = pt[2] - pt[0];
                  break;
                case 1:
                  u = pt[0] - pt[1];
                  v = pt[2] - pt[1];
                  break;
                case 2:
                  u = pt[0] - pt[2];
                  v = pt[1] - pt[2];
                  break;
              }
              return static_cast< OutputVertexNormalComponentType >(
                vcl_acos( u * v ) );
            case AREA:
              return static_cast< OutputVertexNormalComponentType >(
                TriangleType::ComputeArea( pt[0], pt[1], pt[2] ) );
          }
        }
        else
        {
          std::cout <<"Input should be a triangular mesh!!!" <<std::endl;
          return static_cast< OutputVertexNormalComponentType >( 0. );
        }
      }
      else
        return static_cast< OutputVertexNormalComponentType >( 0. );
    }
  }

  template< class TInputMesh, class TOutputMesh >
  void QENormalFilter< TInputMesh, TOutputMesh >::GenerateData( )
  {
    Superclass::GenerateData( );

    ComputeAllFaceNormals( );
    ComputeAllVertexNormals( );
  }
}

#endif
