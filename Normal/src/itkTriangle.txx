#include "itkTriangle.h"

#include <vnl/vnl_math.h>

namespace itk
{
    template< class TPoint >
    Triangle< TPoint >::Triangle( )
    {}

    template< class TPoint >
    Triangle< TPoint >::~Triangle( )
    {}

    template< class TPoint >
    void
    Triangle< TPoint >::SetPoints( const PointType& iP1,
                                   const PointType& iP2,
                                   const PointType& iP3 )
    {
        m_Point[0] = iP1;
        m_Point[1] = iP2;
        m_Point[2] = iP3;
    }

    template< class TPoint >
    void
    Triangle< TPoint >::SetPoint( const unsigned int& iId, const PointType& iPt )
    {
        assert( iId < 3 );
        m_Point[iId] = iPt;
    }

    template< class TPoint >
    typename Triangle< TPoint >::PointType
    Triangle< TPoint >::GetPoint( const unsigned int& iId ) const
    {
        assert( iId < 3 );
        return m_Point[iId];
    }

    template< class TPoint >
    typename Triangle< TPoint >::VectorType
    Triangle< TPoint >::ComputeNormal( ) const
    {
        CrossVectorType cross;
        VectorType w = cross( m_Point[1] - m_Point[0], m_Point[2] - m_Point[0] );
        w.Normalize( );
        return w;
    }

    template< class TPoint >
    typename Triangle< TPoint >::CoordRepType
    Triangle< TPoint >::Cotangent( ) const
    {
        return Triangle< TPoint >::Cotangent( m_Point[0], m_Point[1], m_Point[2] );
    }

    template< class TPoint >
    typename Triangle< TPoint >::PointType
    Triangle< TPoint >::ComputeBarycenter( const CoordRepType& iA,
                                           const CoordRepType& iB,
                                           const CoordRepType& iC ) const
    {
        PointType oPt;
        CoordRepType W = iA + iB + iC;
        if( W != 0. )
        {
            W = 1. / W;

            for( unsigned int dim = 0; dim < PointDimension; dim++ )
                oPt[dim] = W * ( iA * m_Point[0][dim] + iB * m_Point[1][dim] +
                        iC * m_Point[2][dim] );
        }

        return oPt;
    }

    template< class TPoint >
    typename Triangle< TPoint >::PointType
    Triangle< TPoint >::ComputeGravityCenter( ) const
    {
        return ComputeBarycenter( 1., 1., 1. );
    }

    template< class TPoint >
    typename Triangle< TPoint >::PointType
    Triangle< TPoint >::ComputeOrthoCenter( ) const
    {
        CoordRepType a = m_Point[1].SquaredEuclideanDistanceTo( m_Point[2] );
        CoordRepType b = m_Point[0].SquaredEuclideanDistanceTo( m_Point[2] );
        CoordRepType c = m_Point[1].SquaredEuclideanDistanceTo( m_Point[0] );

        CoordRepType Weight[3];
        Weight[0] = ( a + b - c ) * ( c + a - b );
        Weight[1] = ( b + c - a ) * ( a + b - c );
        Weight[2] = ( c + a - b ) * ( b + c - a );

        return ComputeBarycenter( Weight[0], Weight[1], Weight[2] );
    }

    template< class TPoint >
    typename Triangle< TPoint >::PointType
    Triangle< TPoint >::ComputeIncenter( ) const
    {
        CoordRepType a = m_Point[1].EuclideanDistanceTo( m_Point[2] );
        CoordRepType b = m_Point[0].EuclideanDistanceTo( m_Point[2] );
        CoordRepType c = m_Point[1].EuclideanDistanceTo( m_Point[0] );

        return ComputeBarycenter( a, b, c );
    }

    template< class TPoint >
    typename Triangle< TPoint >::CoordRepType
    Triangle< TPoint >::ComputeInRadius( ) const
    {
        CoordRepType a = m_Point[1].EuclideanDistanceTo( m_Point[2] );
        CoordRepType b = m_Point[0].EuclideanDistanceTo( m_Point[2] );
        CoordRepType c = m_Point[1].EuclideanDistanceTo( m_Point[0] );

        CoordRepType t = vnl_math_abs( ( ( b + c - a ) * ( c + a - b ) * ( a + b - c ) ) / ( a + b + c ) );
        return 0.5 * sqrt( t );
    }

    template< class TPoint >
    typename Triangle< TPoint >::CoordRepType
    Triangle< TPoint >::ComputeCircumRadius( ) const
    {
        CoordRepType a = m_Point[1].EuclideanDistanceTo( m_Point[2] );
        CoordRepType b = m_Point[0].EuclideanDistanceTo( m_Point[2] );
        CoordRepType c = m_Point[1].EuclideanDistanceTo( m_Point[0] );
        CoordRepType s = 0.5 * ( a + b + c );
        CoordRepType t = vnl_math_abs( s * ( a + b - s ) * ( a + c - s ) * ( b + c - s ) );

        return ( a * b * c ) / sqrt( t );
    }


    template< class TPoint >
    typename Triangle< TPoint >::PointType
    Triangle< TPoint >::ComputeCircumCenter( ) const
    {
        CoordRepType a = m_Point[1].SquaredEuclideanDistanceTo( m_Point[2] );
        CoordRepType b = m_Point[0].SquaredEuclideanDistanceTo( m_Point[2] );
        CoordRepType c = m_Point[1].SquaredEuclideanDistanceTo( m_Point[0] );

        CoordRepType Weight[3];
        Weight[0] = a * ( b + c - a );
        Weight[1] = b * ( c + a - b );
        Weight[2] = c * ( a + b - c );

        return ComputeBarycenter( Weight[0], Weight[1], Weight[2] );
    }

    template< class TPoint >
    typename Triangle< TPoint >::PointType
    Triangle< TPoint >::ComputeConstrainedCircumCenter( ) const
    {
        PointType oPt;
        CoordRepType a = m_Point[1].SquaredEuclideanDistanceTo( m_Point[2] );
        CoordRepType b = m_Point[0].SquaredEuclideanDistanceTo( m_Point[2] );
        CoordRepType c = m_Point[1].SquaredEuclideanDistanceTo( m_Point[0] );

        CoordRepType Weight[3];
        Weight[0] = a * ( b + c - a );
        Weight[1] = b * ( c + a - b );
        Weight[2] = c * ( a + b - c );

        for( unsigned int i = 0; i < 3; i++ )
        {
            if( Weight[i] < 0. )
                Weight[i] = 0.;
        }

        return ComputeBarycenter( Weight[0], Weight[1], Weight[2] );
    }

    template< class TPoint >
    typename Triangle< TPoint >::CoordRepType
    Triangle< TPoint >::ComputePerimeter( ) const
    {
        CoordRepType a = m_Point[1].EuclideanDistanceTo( m_Point[2] );
        CoordRepType b = m_Point[0].EuclideanDistanceTo( m_Point[2] );
        CoordRepType c = m_Point[1].EuclideanDistanceTo( m_Point[0] );

        return a + b + c;
    }

    template< class TPoint >
    typename Triangle< TPoint >::CoordRepType
    Triangle< TPoint >::ComputeArea( ) const
    {
        CoordRepType a = m_Point[1].EuclideanDistanceTo( m_Point[2] );
        CoordRepType b = m_Point[0].EuclideanDistanceTo( m_Point[2] );
        CoordRepType c = m_Point[1].EuclideanDistanceTo( m_Point[0] );

        CoordRepType s = 0.5 * ( a + b + c );
        return sqrt( s * ( s - a ) * ( s - b ) * ( s - c ) );
    }

    template< class TPoint >
    void
    Triangle< TPoint >::PrintSelf( std::ostream& os, Indent indent ) const
    {
        Superclass::PrintSelf( os, indent );

        os << indent << "Point 1: " << m_Point[0] << std::endl;
        os << indent << "Point 2: " << m_Point[1] << std::endl;
        os << indent << "Point 3: " << m_Point[2] << std::endl;
    }

}
