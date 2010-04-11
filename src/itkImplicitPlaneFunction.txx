#ifndef __itkImplicitPlaneFunction_txx
#define __itkImplicitPlaneFunction_txx

#include "itkImplicitFunctionBase.h"

namespace itk
{

template< typename TOutput, unsigned int VPointDimension, class TInput >
ImplicitPlaneFunction< TOutput, VPointDimension, TInput >::
ImplicitPlaneFunction() : Superclass()
{
  m_Coefficients.Fill( static_cast< OutputType >( 0. ) );
}

template< typename TOutput, unsigned int VPointDimension, class TInput >
ImplicitPlaneFunction< TOutput, VPointDimension, TInput >::
~ImplicitPlaneFunction()
{}

template< typename TOutput, unsigned int VPointDimension, class TInput >
typename ImplicitPlaneFunction< TOutput, VPointDimension, TInput >::OutputType
ImplicitPlaneFunction< TOutput, VPointDimension, TInput >::
Evaluate( const InputType& iPt )
{
  unsigned int k = 0;
  OutputType oValue = static_cast< OutputType >( 0. );

  for( unsigned int dim = 0; dim < PointDimension; ++dim, ++k )
    {
    oValue += m_Coefficients[k] * iPt[dim];
    }
  oValue += m_Coefficients[k];
  
  return oValue;
}

template< typename TOutput, unsigned int VPointDimension, class TInput >
typename ImplicitPlaneFunction< TOutput, VPointDimension, TInput >::GradientType
ImplicitPlaneFunction< TOutput, VPointDimension, TInput >::
Gradient( const InputType& iPt )
{
  GradientType oG;
  for( unsigned int dim = 0; dim < PointDimension; ++dim )
    {
    oG[dim] = m_Coefficients[dim];
    }

  return oG;
}

template< typename TOutput, unsigned int VPointDimension, class TInput >
typename ImplicitPlaneFunction< TOutput, VPointDimension, TInput >::HessianType
ImplicitPlaneFunction< TOutput, VPointDimension, TInput >::
Hessian( const InputType& iPt )
{
  HessianType oH;
  oH.Fill( 0. );

  return oH;
}
}

#endif
