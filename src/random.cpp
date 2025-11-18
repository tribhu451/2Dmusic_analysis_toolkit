#include "random.h"



random_gen::random_gen() : rand_uniform_dist(0.,1.) { 
  seed = ran_dev();
  //seed = 12345 ;  // avoid the randomness by turning on this one (only for check purpose) ... :)
  ran_generator = std::unique_ptr<std::mt19937>(new std::mt19937(seed));

}

random_gen::~random_gen(){
}


double random_gen::rand_uniform(){
  return(rand_uniform_dist(*ran_generator));
}


int random_gen::get_seed() {
  return(seed);
}


// This function generates the random numbers with gaussian 
// distribution having mean=mu and std. dev. = sigma
double random_gen::random_gaussian(double mu, double sigma){
  double x = 0. ; 

  if(mu<0.00000){
    sigma *= -1 ; 
  }
  double limit_of_call = mu + 4.0 * sigma ; 

  int ii=0 ; 
  for( int i=0; i<1000; i++ ){
    ii += 1 ; 
    double test_rand = rand_uniform();
    x = 2.0 * limit_of_call * rand_uniform() - limit_of_call ;
    double gauss = exp( -0.5 * pow( (x - mu) / sigma , 2 ) ) ;  
    if( test_rand < gauss ){
      break ; 
    }
  }

  if(ii > 997 ){
    std::cout << "maximum calls for normal" 
              << "distribution sampling ... " << std::endl ;
    std::cout << "exiting ... " << std::endl ;  
    exit(1);
  }
  
  return x ; 
 
}


