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
  
  PointType p;
  p.Fill( 1. );

  f->Evaluate( p );
  f->Gradient( p );
  f->Hessian( p );

  return EXIT_SUCCESS;
}

