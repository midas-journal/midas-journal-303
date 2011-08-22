INCLUDE_DIRECTORIES( BEFORE
  ${ITK_INCLUDE_DIR}
  "${PROJECT_SOURCE_DIR}/src/" )

SET( NORMAL_EXAMPLE NormalFilter )

ADD_EXECUTABLE( ${NORMAL_EXAMPLE} ${NORMAL_EXAMPLE}.cxx )
TARGET_LINK_LIBRARIES( ${NORMAL_EXAMPLE} ITKCommon ITKQuadEdgeMesh )

FOREACH( normal_type 0 1 2 )
  ADD_TEST( NormalFilterTest${normal_type}
      NormalFilter
      ${DATAROOT}/sphere.vtk ${normal_type} )
ENDFOREACH( normal_type )