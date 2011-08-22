#ifndef __itkTriangle_h
#define __itkTriangle_h

#include <itkDataObject.h>
#include <itkObjectFactory.h>

#include "itkCross.h"

namespace itk
{
  /** \class Triangle
   * \brief Convenient class for various triangles elements computation in
   * 2D or 3D
   * \author Arnaud GELAS
   */
template< typename TPoint >
class Triangle : public DataObject
{
public:
  typedef Triangle Self;
  typedef SmartPointer< Self > Pointer;
  typedef SmartPointer< const Self > ConstPointer;
  typedef DataObject Superclass;

  typedef TPoint PointType;
  typedef typename PointType::CoordRepType CoordRepType;
  typedef typename PointType::VectorType VectorType;
  typedef Cross< VectorType > CrossVectorType;

  itkStaticConstMacro ( PointDimension, unsigned int, PointType::PointDimension );

  itkNewMacro ( Self );

  itkTypeMacro ( Triangle, DataObject );

  void SetPoints ( const PointType& iP1,
                   const PointType& iP2,
                   const PointType& iP3 );

  void SetPoint ( const unsigned int& iId, const PointType& iPt );
  PointType GetPoint ( const unsigned int& iId ) const;

  bool IsObtuse( ) const
  {
    VectorType v01 = m_Point[1] - m_Point[0];
    VectorType v02 = m_Point[2] - m_Point[0];
    VectorType v12 = m_Point[2] - m_Point[1];

    if( v01 * v02 < 0. )
      return true;
    else
      {
      if( v02 * v12 < 0. )
        return true;
      else
        {
        if( v01 * -v12 < 0. )
          return true;
        else
          return false;
        }
      }
  }

  static bool IsObtuse( const PointType& iA,
    const PointType& iB,
    const PointType& iC )
  {
    VectorType v01 = iB - iA;
    VectorType v02 = iC - iA;
    VectorType v12 = iC - iB;

    if( v01 * v02 < 0. )
      return true;
    else
      {
      if( v02 * v12 < 0. )
        return true;
      else
        {
        if( v01 * -v12 < 0. )
          return true;
        else
          return false;
        }
      }
  }
  
  VectorType ComputeNormal( ) const;
  static VectorType ComputeNormal ( const PointType& iA,
    const PointType& iB,
    const PointType& iC )
  {
    CrossVectorType cross;
    VectorType w = cross ( iB - iA, iC - iA );
    w.Normalize( );

    return w;
  }

  CoordRepType Cotangent( ) const;
  static CoordRepType Cotangent ( const PointType& iA,
                                  const PointType& iB,
                                  const PointType& iC )
  {
    VectorType v21 = iA - iB;
    v21.Normalize();
    
    VectorType v23 = iC - iB;
    v23.Normalize();
    
    CoordRepType bound( 0.999999 );

    CoordRepType cos_theta = vnl_math_max( -bound, 
      vnl_math_min( bound, v21 * v23 ) );
            
    return 1. / tan( acos( cos_theta ) );
  }

  PointType ComputeBarycenter ( const CoordRepType& iA,
			        const CoordRepType& iB,
			        const CoordRepType& iC ) const;

  static PointType ComputeBarycenter (
      const CoordRepType& iA1, const PointType& iP1,
      const CoordRepType& iA2, const PointType& iP2,
      const CoordRepType& iA3, const PointType& iP3 )
  {
    PointType oPt;

    for ( unsigned int dim = 0; dim < PointDimension; dim++ )
      oPt[dim] = iA1 * iP1[dim] + iA2 * iP2[dim] + iA3 * iP3[dim];

    return oPt;
  }

  static CoordRepType ComputeAngle( const PointType& iP1, const PointType& iP2,
      const PointType& iP3 )
  {
    VectorType v21 = iP1 - iP2;
    VectorType v23 = iP3 - iP2;

    v21.Normalize();
    v23.Normalize();

    return vcl_acos( v21 * v23 );
  }

  PointType ComputeGravityCenter( ) const;

  static PointType ComputeGravityCenter (
      const PointType& iP1,
      const PointType& iP2,
      const PointType& iP3 )
  {
    PointType oPt;
    CoordRepType inv_3 = 1. / 3.;

    for ( unsigned int dim = 0; dim < PointDimension; dim++ )
      oPt[dim] = ( iP1[dim] + iP2[dim] + iP3[dim] ) * inv_3;

    return oPt;
  }

  PointType ComputeOrthoCenter( ) const;
  PointType ComputeIncenter( ) const;

  CoordRepType ComputeInRadius( ) const;
  CoordRepType ComputeCircumRadius( ) const;

  PointType ComputeCircumCenter( ) const;

  static PointType ComputeCircumCenter (
      const PointType& iP1,
      const PointType& iP2,
      const PointType& iP3 )
  {
    PointType oPt;
    oPt.Fill ( 0. );

    CoordRepType a = iP2.SquaredEuclideanDistanceTo ( iP3 );
    CoordRepType b = iP1.SquaredEuclideanDistanceTo ( iP3 );
    CoordRepType c = iP2.SquaredEuclideanDistanceTo ( iP1 );

    CoordRepType Weight[3];
    Weight[0] = a * ( b + c - a );
    Weight[1] = b * ( c + a - b );
    Weight[2] = c * ( a + b - c );

    CoordRepType SumWeight = Weight[0] + Weight[1] + Weight[2];

    if ( SumWeight != 0. )
      {
      SumWeight = 1. / SumWeight;

      for ( unsigned int dim = 0; dim < PointDimension; dim++ )
        oPt[dim] = ( Weight[0] * iP1[dim] +
            Weight[1] * iP2[dim] + Weight[2] * iP3[dim] ) * SumWeight;
      }


    return oPt;
  }

  PointType ComputeConstrainedCircumCenter( ) const;

  static PointType ComputeConstrainedCircumCenter ( const PointType& iP1,
      const PointType& iP2, const PointType& iP3 )
  {
    PointType oPt;
    CoordRepType a = iP2.SquaredEuclideanDistanceTo ( iP3 );
    CoordRepType b = iP1.SquaredEuclideanDistanceTo ( iP3 );
    CoordRepType c = iP2.SquaredEuclideanDistanceTo ( iP1 );

    CoordRepType Weight[3];
    Weight[0] = a * ( b + c - a );
    Weight[1] = b * ( c + a - b );
    Weight[2] = c * ( a + b - c );

    for ( unsigned int i = 0; i < 3; i++ )
      {
      if ( Weight[i] < 0. )
	Weight[i] = 0.;
      }

    CoordRepType SumWeight = Weight[0] + Weight[1] + Weight[2];

    if ( SumWeight != 0. )
      {
      SumWeight = 1. / SumWeight;

      for ( unsigned int dim = 0; dim < PointDimension; dim++ )
	oPt[dim] = ( Weight[0] * iP1[dim] +
	    Weight[1] * iP2[dim] + Weight[2] * iP3[dim] ) * SumWeight;
      }

    return oPt;
  }

  CoordRepType ComputePerimeter( ) const;
  CoordRepType ComputeArea( ) const;

  static CoordRepType ComputeArea ( const PointType& iP1,
      const PointType& iP2,
      const PointType& iP3 )
  {
    CoordRepType a = iP2.EuclideanDistanceTo ( iP3 );
    CoordRepType b = iP1.EuclideanDistanceTo ( iP3 );
    CoordRepType c = iP2.EuclideanDistanceTo ( iP1 );

    CoordRepType s = 0.5 * ( a + b + c );
    return static_cast< CoordRepType > (
        vcl_sqrt ( s * ( s - a ) * ( s - b ) * ( s - c ) ) );
  }

protected:
  Triangle( );
  virtual ~Triangle( );

  PointType m_Point[3];

  void PrintSelf ( std::ostream& os, Indent indent ) const;


private:
  Triangle ( const Self& );
  void operator = ( const Self& );
};
}

#include "itkTriangle.txx"
#endif
