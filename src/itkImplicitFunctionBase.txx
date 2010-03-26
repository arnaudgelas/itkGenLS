#ifndef __itkImplicitFunctionBase_txx
#define __itkImplicitFunctionBase_txx

namespace itk
{
template< typename TOutput, unsigned int VPointDimension, class TInput >
ImplicitFunctionBase< TOutput, VPointDimension, TInput >::
ImplicitFunctionBase() : Superclass()
{}

template< typename TOutput, unsigned int VPointDimension, class TInput >
ImplicitFunctionBase< TOutput, VPointDimension, TInput >::
~ImplicitFunctionBase()
{}

template< typename TOutput, unsigned int VPointDimension, class TInput >
typename ImplicitFunctionBase< TOutput, VPointDimension, TInput >::GradientType
ImplicitFunctionBase< TOutput, VPointDimension, TInput >::
Normal( const InputType& iPt )
{
  GradientType grad = Gradient( iPt );
  OutputType sq_norm = grad.GetSquaredNorm();
  if( sq_norm > vnl_math::eps )
    {
    sq_norm = 1. / vcl_sqrt( sq_norm );
    return grad * sq_norm;
    }
  else
    {
    return grad;
    }
}

template< typename TOutput, unsigned int VPointDimension, class TInput >
typename ImplicitFunctionBase< TOutput, VPointDimension, TInput >::OutputType
ImplicitFunctionBase< TOutput, VPointDimension, TInput >::Laplacian( const InputType& iPt )
{
  HessianType hess = Hessian( iPt );

  OutputType oValue = 0.;

  for( unsigned int i = 0; i < PointDimension; ++i )
    {
    oValue += hess[i][i];
    }
  return oValue;
}

template< typename TOutput, unsigned int VPointDimension, class TInput >
void
ImplicitFunctionBase< TOutput, VPointDimension, TInput >::
ValueAndGradient( const InputType& iPt,
  OutputType& oValue, GradientType& oGrad )
{
  oValue = Evaluate( iPt );
  oGrad = Gradient( iPt );
}

template< typename TOutput, unsigned int VPointDimension, class TInput >
void
ImplicitFunctionBase< TOutput, VPointDimension, TInput >::
ValueGradientAndHessian( const InputType& iPt,
  OutputType& oValue, GradientType& oGrad, HessianType& oHess ) 
{
  ValueAndGradient( iPt, oValue, oGrad );
  oHess = Hessian( iPt );
}

template< typename TOutput, unsigned int VPointDimension, class TInput >
typename ImplicitFunctionBase< TOutput, VPointDimension, TInput >::OutputType
ImplicitFunctionBase< TOutput, VPointDimension, TInput >::
TaubinErrorValue( const InputType& iPt, const OutputType& iThreshold )
{
  OutputType oValue = Evaluate( iPt );

  GradientType grad = Gradient( iPt );

  if( oValue != 0. )
    {
    // Implicit ImplicitFunctionBase value
    oValue = vnl_math_abs( oValue );

    // Norm of Gradient Vector
    OutputType GradNorm2 = grad.GetSquaredNorm( );

    if( GradNorm2 < iThreshold * iThreshold )
      {
      GradNorm2 = iThreshold * iThreshold;
      }

    if( GradNorm2 > vnl_math::eps * vnl_math::eps )
      {
      oValue /= vcl_sqrt( GradNorm2 );
      }
    }

  return oValue;
}

template< typename TOutput, unsigned int VPointDimension, class TInput >
bool 
ImplicitFunctionBase< TOutput, VPointDimension, TInput >::
IsInside( const InputType& iPt )
{
  OutputType value = Evaluate( iPt );
  return ( value < 0. );
}

template< typename TOutput, unsigned int VPointDimension, class TInput >
bool
ImplicitFunctionBase< TOutput, VPointDimension, TInput >::
IsOutside( const InputType& iPt )
{
  return !this->IsInside( iPt );
}

template< typename TOutput, unsigned int VPointDimension, class TInput >
typename ImplicitFunctionBase< TOutput, VPointDimension, TInput >::OutputType
ImplicitFunctionBase< TOutput, VPointDimension, TInput >::
MeanCurvature( const InputType& iPt )
{
  return MeanCurvature( Dispatch< PointDimension >(), iPt );
}

template< typename TOutput, unsigned int VPointDimension, class TInput >
typename ImplicitFunctionBase< TOutput, VPointDimension, TInput >::OutputType
ImplicitFunctionBase< TOutput, VPointDimension, TInput >::
GaussKroneckerCurvature( const InputType& iPt )
{
  return GaussKroneckerCurvature( Dispatch< PointDimension >(), iPt );
}

template< typename TOutput, unsigned int VPointDimension, class TInput >
typename ImplicitFunctionBase< TOutput, VPointDimension, TInput >::OutputType
ImplicitFunctionBase< TOutput, VPointDimension, TInput >::
MeanCurvature( const Dispatch< 3 >&,
  const InputType& iPt )
{
  GradientType grad = Gradient( iPt );
  HessianType hessian = Hessian( iPt );

  OutputType norm_g = grad.GetNorm();

  if( norm_g < vnl_math::eps )
    {
    return 0.;
    }
  else
    {
    HessianType Identity;
    Identity.SetIdentity();
    OutputType inv_norm_g = 1. / norm_g;
    grad *= inv_norm_g;
    OutputType factor = ( Identity - grad * grad.transpose() ) * hessian;

    return -inv_norm_g * factor;
    }
}

template< typename TOutput, unsigned int VPointDimension, class TInput >
typename ImplicitFunctionBase< TOutput, VPointDimension, TInput >::OutputType
ImplicitFunctionBase< TOutput, VPointDimension, TInput >::
MeanCurvature( const DispatchBase&,
  const InputType& iPt )
{
  GradientType grad = Gradient( iPt );
  HessianType  hessian = Hessian( iPt );

  OutputType sq_norm_g = grad.GetSquaredNorm();
  OutputType norm_g = vcl_sqrt( sq_norm_g );
  OutputType gHgt = grad * hessian * grad.transpose();
  OutputType trace_H = 0.;

  for( unsigned int dim = 0; dim < PointDimension; ++dim )
    {
    trace_H += hessian[dim][dim];
    }

  return ( 1. / ( ( PointDimension - 1 ) * norm_g * sq_norm_g ) *
    (gHgt - sq_norm_g * trace_H) );
}

template< typename TOutput, unsigned int VPointDimension, class TInput >
typename ImplicitFunctionBase< TOutput, VPointDimension, TInput >::OutputType
ImplicitFunctionBase< TOutput, VPointDimension, TInput >::
GaussKroneckerCurvature( const Dispatch< 2 >&,
  const InputType& iPt )
{
  return 0.;
}

template< typename TOutput, unsigned int VPointDimension, class TInput >
typename ImplicitFunctionBase< TOutput, VPointDimension, TInput >::OutputType
ImplicitFunctionBase< TOutput, VPointDimension, TInput >::
GaussKroneckerCurvature( const DispatchBase&,
  const InputType& iPt )
{
  GradientType g = Gradient( iPt );
  OutputType sq_norm_g = g.GetSquaredNorm();
  if( sq_norm_g == 0. )
    {
    return 0.;
    }
  else
    {
    OutputType norm_g = vcl_sqrt( sq_norm_g );
    HessianType H = Hessian( iPt );
    OutputType coeff = ( PointDimension % 2 == 1 ) ? -norm_g : norm_g;

    for( unsigned int dim = 0; dim < PointDimension; ++dim )
      {
      coeff *= norm_g;
      }
  // should be conjuguate but H is real...
  OutputType gtH_g = g.transpose() * H.transpose() * g;

  return gtH_g / coeff;
  }
}

template< typename TOutput, unsigned int VPointDimension, class TInput >
typename ImplicitFunctionBase< TOutput, VPointDimension, TInput >::HessianType
ImplicitFunctionBase< TOutput, VPointDimension, TInput >::
CurvatureTensor( const Dispatch<3>&, const InputType& iPt )
{
  GradientType g = Gradient( iPt );
  HessianType H = Hessian( iPt );

  OutputType inv_sq_norm = 1. / g.GetSquaredNorm();
  OutputType inv_norm = vcl_sqrt( inv_sq_norm );

  Matrix< OutputType, PointDimension, 1 > Gradient;

  for( unsigned int i = 0; i < PointDimension; i++ )
    {
    Gradient[i] = g[i];
    }

  ImplicitFunctionBaseGradient2( iPt, Hessian );

  HessianType oTensor = ( H -
    inv_sq_norm * g.transpose() * H * g ) * inv_norm;

  return oTensor;
}

template< typename TOutput, unsigned int VPointDimension, class TInput >
typename ImplicitFunctionBase< TOutput, VPointDimension, TInput >::HessianType
ImplicitFunctionBase< TOutput, VPointDimension, TInput >::
CurvatureTensor( const DispatchBase&, const InputType& iPt )
{
  HessianType oTensor;
  oTensor.Fill( 0. );
  return oTensor;
}

template< typename TOutput, unsigned int VPointDimension, class TInput >
void
ImplicitFunctionBase< TOutput, VPointDimension, TInput >::
PrincipalCurvatures( const DispatchBase&,
  const InputType& iPt,
  OutputType& oKmin,
  OutputType& oKmax,
  GradientType& oEmin,
  GradientType& oEmax )
{
  oKmin = 0.;
  oKmax = 0.;
  oEmin.Fill( 0. );
  oEmax.Fill( 0. );
}

template< typename TOutput, unsigned int VPointDimension, class TInput >
void
ImplicitFunctionBase< TOutput, VPointDimension, TInput >::
PrincipalCurvatures( const Dispatch<3>&,
  const InputType& iPt,
  OutputType& oKmin,
  OutputType& oKmax,
  GradientType& oEmin,
  GradientType& oEmax )
{
  OutputType mean_curvature = MeanCurvature( Dispatch<3>(), iPt );
  OutputType gauss_curvature =
    GaussKroneckerCurvature( Dispatch<3>(), iPt );

  OutputType delta =
    vnl_math_abs( mean_curvature * mean_curvature - gauss_curvature );

  OutputType sqrt_delta = vcl_sqrt( delta );

  oKmax = mean_curvature + sqrt_delta;
  oKmin = mean_curvature - sqrt_delta;

  if( sqrt_delta < 1e-6 )
    {
    oEmin.Fill( 0. );
    oEmax.Fill( 0. );
    }
  else
    {
    HessianType tensor = CurvatureTensor( Dispatch<3>(), iPt );

    /* solve for eigenvectors */
    OutputType t[3];
    t[0] = tensor[0][0] + oKmax;
    t[1] = tensor[1][1] + oKmax;
    t[2] = tensor[2][2] + oKmax;

    GradientType v[3];
    OutputType len[3];

    v[0][0] = tensor[0][1] * tensor[1][2] - tensor[0][2] * t[1];
    v[1][0] = tensor[0][2] * tensor[1][0] - tensor[1][2] * t[0];
    v[2][0] = t[0] * t[1] - tensor[0][1] * tensor[1][0];
    len[0] = vcl_sqrt( v[0][0] * v[0][0] +
      v[1][0] * v[1][0] + v[2][0] * v[2][0] );

    v[0][1] = tensor[0][1] * t[2] - tensor[0][2] * tensor[2][1];
    v[1][1] = tensor[0][2] * tensor[2][0]- t[0] * t[2];
    v[2][1] = t[0] * tensor[2][1] - tensor[0][1] * tensor[2][0];
    len[1] = vcl_sqrt( v[0][1] * v[0][1] +
      v[1][1] * v[1][1] + v[2][1] * v[2][1] );

    v[0][2] = t[1] * t[2] - tensor[1][2] * tensor[2][1];
    v[1][2] = tensor[1][2] * tensor[2][0]- tensor[1][0] * t[2];
    v[2][2] = tensor[1][0] * tensor[2][1] - tensor[2][0]* t[1];
    len[2] = vcl_sqrt( v[0][2] * v[0][2] +
      v[1][2] * v[1][2] + v[2][2] * v[2][2] );

    unsigned int index = 0;
    OutputType len_max = len[index];

    if( len[1] > len_max )
      {
      index = 1;
      len_max = len[1];
      }
    if( len[2] > len_max )
      {
      index = 2;
      len_max = len[2];
      }

    oEmax = v[index] / len[index];

    GradientType n = Normal( iPt );

    oEmin[0] = oEmax[1] * n[2] - oEmax[2] * n[1];
    oEmin[1] = oEmax[2] * n[0] - oEmax[0] * n[2];
    oEmin[2] = oEmax[0] * n[1] - oEmax[1] * n[0];
    }
}

template< typename TOutput, unsigned int VPointDimension, class TInput >
void
ImplicitFunctionBase< TOutput, VPointDimension, TInput >::
PrincipalCurvaturesTensor( const InputType& iPt,
  OutputType& oKmin,
  OutputType& oKmax,
  GradientType& oEmin,
  GradientType& oEmax )
{
  PrincipalCurvaturesTensor( Dispatch< PointDimension >(), iPt,
    oKmin, oKmax, oEmin, oEmax );
}

template< typename TOutput, unsigned int VPointDimension, class TInput >
typename ImplicitFunctionBase< TOutput, VPointDimension, TInput >::HessianType
ImplicitFunctionBase< TOutput, VPointDimension, TInput >::CurvatureTensor( const InputType& iPt )
{
  return CurvatureTensor( Dispatch< PointDimension >(), iPt );
}

template< typename TOutput, unsigned int VPointDimension, class TInput >
void
ImplicitFunctionBase< TOutput, VPointDimension, TInput >::
PrintSelf(std::ostream& os, Indent indent) const
{
  Superclass::PrintSelf( os, indent );
}

}
#endif
