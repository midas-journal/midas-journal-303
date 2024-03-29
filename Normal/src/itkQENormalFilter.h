#ifndef __itkQENormalFilter_h
#define __itkQENormalFilter_h

#include <itkQuadEdgeMeshToQuadEdgeMeshFilter.h>
#include <itkQuadEdgeMeshPolygonCell.h>
#include "itkTriangle.h"

namespace itk
  {
  /** \brief Filter which computes normals to faces and vertices and store it in
  * the output mesh.
  * \todo Fix run-time issues regarding the difference between the Traits of
  * TInputMesh and the one of TOutputMesh. Right now, it only works if
  * TInputMesh::MeshTraits == TOutputMesh::MeshTraits
  * (and of course it requires that the output have some itk::Vector for point
  * data and cell data.
  */
  template< class TInputMesh, class TOutputMesh >
  class QENormalFilter :
    public QuadEdgeMeshToQuadEdgeMeshFilter< TInputMesh, TOutputMesh >
  {
    public:
      typedef QENormalFilter Self;
      typedef QuadEdgeMeshToQuadEdgeMeshFilter< TInputMesh, TOutputMesh >
      Superclass;
      typedef SmartPointer< Self > Pointer;
      typedef SmartPointer< const Self > ConstPointer;

      itkNewMacro ( Self );

      itkTypeMacro ( QENormalFilter,
                     QuadEdgeMeshToQuadEdgeMeshFilter );

      typedef TInputMesh InputMeshType;
      typedef typename InputMeshType::Pointer InputMeshPointer;
      typedef typename InputMeshType::PointIdentifier InputPointIdentifier;
      typedef typename InputMeshType::PointType InputPointType;
      typedef typename InputMeshType::VectorType InputVectorType;
      typedef typename InputMeshType::QEType InputQEType;

      typedef TOutputMesh OutputMeshType;
      typedef typename OutputMeshType::Pointer OutputMeshPointer;
      typedef typename OutputMeshType::PointType OutputPointType;
      typedef typename OutputPointType::VectorType OutputVectorType;
      typedef typename OutputMeshType::QEType OutputQEType;
      typedef typename OutputMeshType::PointIdentifier OutputPointIdentifier;
      typedef typename OutputMeshType::PointIdIterator OutputPointIdIterator;
      typedef typename OutputMeshType::PointsContainerPointer
        OutputPointsContainerPointer;
      typedef typename OutputMeshType::PointsContainerIterator
        OutputPointsContainerIterator;
      typedef typename OutputMeshType::CellType OutputCellType;
      typedef typename OutputMeshType::CellIdentifier OutputCellIdentifier;
      typedef typename OutputMeshType::CellAutoPointer OutputCellAutoPointer;
      typedef typename OutputMeshType::CellsContainerConstIterator
        OutputCellsContainerPointer;
      typedef typename OutputMeshType::CellsContainerConstIterator
        OutputCellsContainerConstIterator;

      typedef Triangle< OutputPointType > TriangleType;

      typedef QuadEdgeMeshPolygonCell< OutputCellType > OutputPolygonType;
      typedef typename OutputPolygonType::SelfAutoPointer
        OutputPolygonAutoPointer;

      typedef typename OutputMeshType::CellDataContainer OutputCellDataContainer;
      typedef typename OutputMeshType::PointDataContainer
        OutputPointDataContainer;

      typedef typename OutputMeshType::MeshTraits OutputMeshTraits;
      typedef typename OutputMeshTraits::PixelType OutputVertexNormalType;
      typedef typename OutputVertexNormalType::ValueType
        OutputVertexNormalComponentType;

      typedef typename OutputMeshTraits::CellPixelType OutputFaceNormalType;
      typedef typename OutputFaceNormalType::ValueType
        OutputFaceNormalComponentType;

      enum WeightType
      {
        GOURAUD = 0, // Uniform weights
        THURMER, // Angle on a triangle at the given vertex
        AREA
      };

      itkSetMacro ( Weight, WeightType );
      itkGetMacro ( Weight, WeightType );

    protected:
      QENormalFilter( );
      ~QENormalFilter( );

      WeightType m_Weight;

      /** \brief Compute the normal to a face iPoly. It assumes that iPoly != 0
      * and
      * iPoly is a Triangle, i.e. 3 points only.
      * \note The normal computation itself can be further improved by making
      * possible to cast a CellType into a TriangleType.
      */
      OutputFaceNormalType ComputeFaceNormal ( OutputPolygonType* iPoly );

      /** \brief Compute the normal to all faces on the mesh.
      * \note This method should be implemented in a multi-thread way in order
      * to reduce the processing time.
      */
      void ComputeAllFaceNormals( );

      /** \brief Compute the normal to all vertices on the mesh.
      * \note This method should be implemented in a multi-thread way in order
      * to reduce the processing time.
      */
      void ComputeAllVertexNormals( );

      /** \brief Compute the normal to one vertex by a weighted sum of the faces
      * normal in the 0-ring.
      * \note The weight is chosen by the member m_Weight.
      */
      OutputVertexNormalType ComputeVertexNormal (
        const OutputPointIdentifier& iId );

      /** \brief Definition of the weight in the 0-ring used for the vertex
      * normal computation. By default m_Weight = THURMER;
      */
      OutputVertexNormalComponentType Weight ( const OutputPointIdentifier& iPId,
          const OutputCellIdentifier& iCId );

      /** \note Calling Superclass::GenerateData( ) is the longest part in the
      * filter! Something must be done in the class
      * itkQuadEdgeMeshToQuadEdgeMeshFilter.
      */
      void GenerateData( );

    private:
      QENormalFilter ( const Self& );
      void operator = ( const Self& );
    };
}

#include "itkQENormalFilter.txx"
#endif
