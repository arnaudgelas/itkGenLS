#ifndef __itkImplicitQuadraticFunction_txx
#define __itkImplicitQuadraticFunction_txx

#include "itkImplicitFunctionBase.h"

namespace itk
{

template< typename TOutput, unsigned int VPointDimension, class TInput >
ImplicitQuadraticFunction< TOutput, VPointDimension, TInput >::
ImplicitQuadraticFunction() : Superclass()
{
  m_Coefficients.Fill( static_cast< OutputType >( 0. ) );
}

template< typename TOutput, unsigned int VPointDimension, class TInput >
ImplicitQuadraticFunction< TOutput, VPointDimension, TInput >::
~ImplicitQuadraticFunction()
{}

template< typename TOutput, unsigned int VPointDimension, class TInput >
typename ImplicitQuadraticFunction< TOutput, VPointDimension, TInput >::OutputType
ImplicitQuadraticFunction< TOutput, VPointDimension, TInput >::
Evaluate( const InputType& iPt )
{
  OutputType oValue = static_cast< OutputType >( 0. );
  
  unsigned int dim, dim2, k = 0;
  
  for( dim = 0; dim < PointDimension; ++dim )
    {
    for( dim2 = dim; dim2 < PointDimension; ++dim2 )
      {
      oValue += m_Coefficients[k++] * 
        static_cast< OutputType >( iPt[dim] ) * 
        static_cast< OutputType >( iPt[dim2] );
      }
    }

  for( dim = 0; dim < PointDimension; ++dim )
    {
    oValue += m_Coefficients[k++] * static_cast< OutputType >( iPt[dim] );
    }

  oValue += m_Coefficients[k++];
  
  return oValue;
}

template< typename TOutput, unsigned int VPointDimension, class TInput >
typename ImplicitQuadraticFunction< TOutput, VPointDimension, TInput >::GradientType
ImplicitQuadraticFunction< TOutput, VPointDimension, TInput >::
Gradient( const InputType& iPt )
{
  GradientType oG;
  oG.Fill( 0. );

  unsigned int dim, dim2, k = 0;

  for( dim = 0; dim < PointDimension; ++dim )
    {
    oG[dim] += static_cast< OutputType >( 2.) * m_Coefficients[k++] * 
      static_cast< OutputType >( iPt[dim] );

    for( dim2 = dim+1; dim2 < PointDimension; ++dim2 )
      {
      oG[dim] += m_Coefficients[k++] * 
        static_cast< OutputType >( iPt[dim2] );
      }
    }

  for( dim = 0; dim < PointDimension; ++dim )
    {
    oG[dim] += m_Coefficients[k++];
    }

  return oG;
}

template< typename TOutput, unsigned int VPointDimension, class TInput >
typename ImplicitQuadraticFunction< TOutput, VPointDimension, TInput >::HessianType
ImplicitQuadraticFunction< TOutput, VPointDimension, TInput >::
Hessian( const InputType& iPt )
{
  HessianType oH;

  unsigned int dim, dim2, k = 0;

  for( dim = 0; dim < PointDimension; ++dim )
    {
    oH[dim][dim] = static_cast< OutputType >( 2. ) * m_Coefficients[k++];

    for( dim2 = dim+1; dim2 < PointDimension; ++dim2 )
      {
      oH[dim][dim2] = m_Coefficients[k++];
      oH[dim2][dim] = oH[dim][dim2];
      }
    }
  return oH;
}
}

#endif
