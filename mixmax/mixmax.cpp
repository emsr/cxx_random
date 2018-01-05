/*
 * mixmax.cpp version 0.1
 * 
 *
 */

#include <iostream>
#include <exception>

#include "mixmax.hpp"

#define MOD_PAYNE(k) ((((k)) & MERSBASE) + (((k)) >> BITS) )
#define MOD_MERSENNE(k) MOD_PAYNE(k)

myuint mixmax_engine::MOD_MULSPEC(myuint k){
    switch (N) {
        case 17:
            return 0;
            break;
        case 256:
            if(SPECIAL==-1){
                return  (MERSBASE - (k));
            }else{
               return fmodmulM61( 0, SPECIAL , (k) ) ;
            }
            break;
        case 240:
            return fmodmulM61( 0, SPECIAL , (k) );
            break;
        default:
            std::cerr << "MIXMAX ERROR: " << "Disallowed value of parameter N\n";
            break;
    }
}

mixmax_engine::mixmax_engine()
// constructor, with no params, fast and seeds with a unit vector
{
    seed_vielbein(&S,0);
}

mixmax_engine::mixmax_engine(myID_t clusterID, myID_t machineID, myID_t runID, myID_t  streamID)
// constructor, no need to allocate, just seed
{
    seed_uniquestream( &S, clusterID,  machineID,  runID,  streamID );
}

#define MULWU(k) (( (k)<<(SPECIALMUL) & M61) | ( (k) >> (BITS-SPECIALMUL))  )

myuint mixmax_engine::iterate_raw_vec(myuint* Y, myuint sumtotOld){
    // operates with a raw vector, uses known sum of elements of Y
    int i;

    myuint temp2 = Y[1];


    myuint  tempP, tempV;
    Y[0] = ( tempV = sumtotOld);
    myuint sumtot = Y[0], ovflow = 0; // will keep a running sum of all new elements
    tempP = 0;              // will keep a partial sum of all old elements
    for (i=1; i<N; i++){
        if (SPECIALMUL!=0){
        myuint tempPO = MULWU(tempP);
        tempP = modadd(tempP,Y[i]);
        tempV = MOD_MERSENNE(tempV + tempP + tempPO); // edge cases ?
        }else{
        tempP = modadd(tempP , Y[i]);
        tempV = modadd(tempV , tempP);
        }
        Y[i] = tempV;
        sumtot += tempV; if (sumtot < tempV) {ovflow++;}
    }
    if (SPECIAL!=0){
    temp2 = MOD_MULSPEC(temp2);
    Y[2] = modadd( Y[2] , temp2 );
    sumtot += temp2; if (sumtot < temp2) {ovflow++;}
    }
    return MOD_MERSENNE(MOD_MERSENNE(sumtot) + (ovflow <<3 ));
}


//rng_state_t* mixmax_engine::rng_alloc()
//{
//    /* allocate the state */
//    rng_state_t  *p = (rng_state_t*)malloc(sizeof(rng_state_t));
//    p->fh=NULL; // by default, set the output file handle to stdout
//    return p;
//}



myuint mixmax_engine::get_next() {
    int i;
    i=S.counter;
    
    if (i<=(N-1) ){
        S.counter++;
        return S.V[i];
    }else{
        
        if(N==256 && SPECIAL==-1){
            S.sumtot = iterate_raw_vec(S.V.data(), S.sumtot);
            S.sumtot = iterate_raw_vec(S.V.data(), S.sumtot);
            S.sumtot = iterate_raw_vec(S.V.data(), S.sumtot);
            // skipping by 2: Requirement from Lorenzo and Fred
        }else{
            S.sumtot = iterate_raw_vec(S.V.data(), S.sumtot);
        }
        S.counter=2;
        return S.V[1];
    }
}

double mixmax_engine::get_next_float()				// Returns a random double with all 53 bits random, in the range (0,1]
{    /* cast to signed int trick suggested by Andrzej Görlich     */
    int64_t Z=(int64_t)get_next();
    double F;
#if defined(__GNUC__) && (__GNUC__ < 5) && (!defined(__ICC)) && defined(__x86_64__) && defined(__SSE2_MATH__) && defined(USE_INLINE_ASM)
    //#warning Using the inline assembler
    /* using SSE inline assemly to zero the xmm register, just before int64 -> double conversion,
     not necessary in GCC-5 or better, but huge penalty on earlier compilers
     */
    __asm__  __volatile__("pxor %0, %0; "
                          :"=x"(F)
                          );
#endif
    F=Z;
    return F*INV_MERSBASE;
    
}

void mixmax_engine::seed_vielbein(rng_state_t* X, unsigned int index)
{
    int i;
    if (index<N){
        for (i=0; i < N; i++){
            S.V[i] = 0;
        }
        S.V[index] = 1;
    }else{
        //fprintf(stderr, "Out of bounds index, is not ( 0 <= index < N  )\n");
        std::cerr << "MIXMAX ERROR: " << ARRAY_INDEX_OUT_OF_BOUNDS << "Out of bounds index, is not ( 0 <= index < N  )\n";
        std::terminate();
    }
    S.counter = N;  // set the counter to N if iteration should happen right away
    S.sumtot = 1;   //(index ? 1:0);
}


void mixmax_engine::seed_uniquestream( rng_state_t* Xin, myID_t clusterID, myID_t machineID, myID_t runID, myID_t  streamID ){
    seed_vielbein(Xin,0);
    Xin->sumtot = apply_bigskip(Xin->V.data(), Xin->V.data(),  clusterID,  machineID,  runID,   streamID );
//   if (Xin->fh==NULL){Xin->fh=stdout;} // if the filehandle is not yet set, make it stdout
    Xin->counter = 1;
}


myuint mixmax_engine::apply_bigskip( myuint* Vout, myuint* Vin, myID_t clusterID, myID_t machineID, myID_t runID, myID_t  streamID ){
    /*
     makes a derived state vector, Vout, from the mother state vector Vin
     by skipping a large number of steps, determined by the given seeding ID's
     
     it is mathematically guaranteed that the substreams derived in this way from the SAME (!!!) Vin will not collide provided
     1) at least one bit of ID is different
     2) less than 10^100 numbers are drawn from the stream
     (this is good enough : a single CPU will not exceed this in the lifetime of the universe, 10^19 sec,
     even if it had a clock cycle of Planch time, 10^44 Hz )
     
     Caution: never apply this to a derived vector, just choose some mother vector Vin, for example the unit vector by seed_vielbein(X,0),
     and use it in all your runs, just change runID to get completely nonoverlapping streams of random numbers on a different day.
     
     clusterID and machineID are provided for the benefit of large organizations who wish to ensure that a simulation
     which is running in parallel on a large number of  clusters and machines will have non-colliding source of random numbers.
     
     did i repeat it enough times? the non-collision guarantee is absolute, not probabilistic
     
     */
    
    
const	myuint skipMat256old[128][256] =
#include "mixmax_skip_N256.oldS.c"
;
const	myuint skipMat256[128][256] =
#include "mixmax_skip_N256.c"
;
const	myuint skipMat240[128][240] =
#include "mixmax_skip_N240.c"
;
const	myuint skipMat17[128][17] =
#include "mixmax_skip_N17.c"
;

const myuint* skipMat[128];
    switch (N) {
        case 17:
            for (int i=0; i<128; i++) { skipMat[i] = skipMat17[i];}
            break;
        case 240:
            for (int i=0; i<128; i++) { skipMat[i] = skipMat240[i];}
            break;
        case 256:
            if(SPECIAL==-1){
                for (int i=0; i<128; i++) { skipMat[i] = skipMat256old[i];}
            }else{
                for (int i=0; i<128; i++) { skipMat[i] = skipMat256[i];}
            }
            break;
        default:
            exit(-1);
            break;
    }
    
    myID_t IDvec[4] = {streamID, runID, machineID, clusterID};
    int r,i,j,  IDindex;
    myID_t id;
    myuint Y[N], cum[N];
    myuint coeff;
    myuint* rowPtr;
    myuint sumtot=0;
    
    
    for (i=0; i<N; i++) { Y[i] = Vin[i]; sumtot = modadd( sumtot, Vin[i]); } ;
    for (IDindex=0; IDindex<4; IDindex++) { // go from lower order to higher order ID
        id=IDvec[IDindex];
        //printf("now doing ID at level %d, with ID = %d\n", IDindex, id);
        r = 0;
        while (id){
            if (id & 1) {
                rowPtr = (myuint*)skipMat[r + IDindex*8*sizeof(myID_t)];
                //printf("free coeff for row %d is %llu\n", r, rowPtr[0]);
                for (i=0; i<N; i++){ cum[i] = 0; }
                for (j=0; j<N; j++){              // j is lag, enumerates terms of the poly
                    // for zero lag Y is already given
                    coeff = rowPtr[j]; // same coeff for all i
                    for (i =0; i<N; i++){
                        cum[i] =  fmodmulM61( cum[i], coeff ,  Y[i] ) ;
                    }
                    sumtot = iterate_raw_vec(Y, sumtot);
                }
                sumtot=0;
                for (i=0; i<N; i++){ Y[i] = cum[i]; sumtot = modadd( sumtot, cum[i]); } ;
            }
            id = (id >> 1); r++; // bring up the r-th bit in the ID		
        }		
    }
    sumtot=0;
    for (i=0; i<N; i++){ Vout[i] = Y[i]; sumtot = modadd( sumtot, Y[i]); } ;  // returns sumtot, and copy the vector over to Vout
    return (sumtot) ;
}

#if defined(__x86_64__)
inline myuint mixmax_engine::mod128(__uint128_t s){
    myuint s1;
    s1 = ( (  ((myuint)s)&MERSBASE )    + (  ((myuint)(s>>64)) * 8 )  + ( ((myuint)s) >>BITS) );
    return	MOD_MERSENNE(s1);
}

inline myuint mixmax_engine::fmodmulM61(myuint cum, myuint a, myuint b){
    __uint128_t temp;
    temp = (__uint128_t)a*(__uint128_t)b + cum;
    return mod128(temp);
}

#else // on all other platforms, including 32-bit linux, PPC and PPC64, ARM and all Windows
#define MASK32 0xFFFFFFFFULL

inline myuint mixmax_engine::fmodmulM61(myuint cum, myuint s, myuint a)
{
    register myuint o,ph,pl,ah,al;
    o=(s)*a;
    ph = ((s)>>32);
    pl = (s) & MASK32;
    ah = a>>32;
    al = a & MASK32;
    o = (o & M61) + ((ph*ah)<<3) + ((ah*pl+al*ph + ((al*pl)>>32))>>29) ;
    o += cum;
    o = (o & M61) + ((o>>61));
    return o;
}
#endif

myuint mixmax_engine::modadd(myuint foo, myuint bar){
#if (defined(__x86_64__) || defined(__i386__)) &&  defined(__GNUC__) && defined(USE_INLINE_ASM)
    //#warning Using assembler routine in modadd
    myuint out;
    /* Assembler trick suggested by Andrzej Görlich     */
    __asm__ ("addq %2, %0; "
             "btrq $61, %0; "
             "adcq $0, %0; "
             :"=r"(out)
             :"0"(foo), "r"(bar)
             );
    return out;
#else
    return MOD_MERSENNE(foo+bar);
#endif
}

void mixmax_engine::print_state(){ // (std::ostream& ost){
    int j;
    fprintf(stdout, "mixmax state, file version 1.0\n" );
    fprintf(stdout, "N=%u; V[N]={", rng_get_N() );
    for (j=0; (j< (rng_get_N()-1) ); j++) {
        fprintf(stdout, "%llu, ", S.V[j] );
    }
    fprintf(stdout, "%llu", S.V[rng_get_N()-1] );
    fprintf(stdout, "}; " );
    fprintf(stdout, "counter=%u; ", S.counter );
    fprintf(stdout, "sumtot=%llu;\n", S.sumtot );
}

mixmax_engine mixmax_engine::Branch(){
    // suggested by Lorenzo
    mixmax_engine tmp=*this;
    this->BranchMother();
    tmp.BranchDaughter(5);
    tmp.print_state();
    this->print_state();
    return tmp;
}

mixmax_engine& mixmax_engine::operator=(const mixmax_engine& other ){
    S = other.S;
   // S.V = S.A.data();
    return *this;
}

void mixmax_engine::BranchDaughter(int b){ // valid values are between b=5 and b=60
    if(b>60) {std::cerr << "MIXMAX ERROR: " << "Disallowed value of parameter b in BranchDaughter\n"; exit(-1);}
    // Dont forget to branch mother, when you branch the daughter, or else you will have collisions!
    printf("V[N] = %llu, %d\n", S.V[N], b); S.V[N] ^= (1<<(BITS-b)); S.V[N] |= 1; printf("V[N] = %llu\n", S.V[N]);
    S.sumtot = iterate_raw_vec(S.V.data(), S.sumtot); printf("iterating!\n");
}

void mixmax_engine::BranchMother(){
    BranchDaughter(4); // same thing, but b must be different
}