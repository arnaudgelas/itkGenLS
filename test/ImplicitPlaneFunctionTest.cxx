#include "itkPoint.h"
#include "itkImplicitPlaneFunction.h"

int main( int argc, char** argv )
{
  const unsigned int Dimension = 2;
  typedef itk::Point< double, Dimension > PointType;
  typedef float OutputType;

  typedef itk::ImplicitPlaneFunction< OutputType, Dimension, PointType >
    PlaneType;
  PlaneType::Pointer f = PlaneType::New();
 
  PlaneType::CoefficientVectorType coeff;
  coeff.Fill( 1. );
  f->SetCoefficients( coeff );
   
  PointType p;
  p.Fill( 1. );

  if( f->Evaluate( p ) != 3 )
    {
    std::cerr <<"Evaluate(p) Failure" <<std::endl;
    return EXIT_FAILURE;
    }

  PlaneType::GradientType g = f->Gradient( p );
  if( ( g[0] != 1 ) || ( g[1] != 1 ) )
    {
    std::cerr <<"Gradient( p ) Failure" <<std::endl;
    return EXIT_FAILURE;
    }

  PlaneType::HessianType h = f->Hessian( p );
  if( ( h[0][0] != 0 ) || ( h[0][1] != 0 ) ||
      ( h[1][0] != 0 ) || ( h[1][1] != 0 ) )
    {
    std::cerr <<"Hessian( p ) Failure" <<std::endl;
    return EXIT_FAILURE;
    }

  return EXIT_SUCCESS;
}

