import cern.jet.random.Exponential;
import cern.jet.random.Normal;
import cern.jet.random.engine.MersenneTwister;
import cern.jet.random.engine.RandomEngine;

/**
 * Created with IntelliJ IDEA.
 * User: Jayna
 * Date: 17/05/2012
 * Time: 21:59
 * To change this template use File | Settings | File Templates.
 */
public class EpiParams {

    //default parameter values?
    static int N = (int)3e6;
    static double R0 = 2;//3.0;
    static int lifeSpan = 30;  //in years
    static int durationOfInfection = 4; // in days
    static int waningImmunity = 4; //in years

    static double psi = 1.0;
    static double psi_i = 0.5;
    static double psi_j = 0.5;
    static double epsilon = 0.0;
    static double simulationTime = 60;
    static double simulationStartTime = 0.0;
    static double simulationEndTime = simulationTime*365.25;

    static int viralLoad = 1000;
    static int mutnRate = 500; //500 gives 6e-3; 10 gives 4e-4

    static double gamma = 1.0/((double)waningImmunity*365.25);
    static double mig = 0.0;
    static double nu_s = 1.0/((double)durationOfInfection);
    static double nu_co = nu_s;
    static double mu = 1.0/((double)lifeSpan*365.25);
    static double beta_s = R0/((double)durationOfInfection);
    static double beta_co = beta_s;

    //initial conditions - start at equilibrium values:
    static int S_init = (int)Math.floor(N/R0);
    static int Is_init = 490*2;//(int)Math.floor((mu*N)/(beta_s*(R0-1.0)));
    static int Ico_init = 10*2;//(int)Math.floor(beta_co/((N*(nu_co+mu))*(Math.pow(Is_init,2))));
    static int R_init = N-S_init-Is_init-Ico_init;

    static double epsilon_endemic = 0.0;
    static double epsilon_seasonal = 0.15;
    static double m_ij = 0.0005;
    static int Is_init_i = 8000;//(int)Math.floor((mu*N)/(beta_s*(R0-1.0)));;
    static int Is_init_j = 8000;//(int)Math.floor((mu*N)/(beta_s*(R0-1.0)));;
    static int S_init_i = S_init;
    static int S_init_j = S_init;
    static int Ico_init_i = 10; //(int)Math.floor(beta_co/((N*(nu_co+mu))*(Math.pow(Is_init_i,2))));;
    static int Ico_init_j = 10;//(int)Math.floor(beta_co/((N*(nu_co+mu))*(Math.pow(Is_init_j,2))));;
    static int R_init_i = N-S_init_i-Is_init_i-Ico_init_i;
    static int R_init_j = N-S_init_j-Is_init_j-Ico_init_j;

    static RandomEngine randomGenerator = new MersenneTwister();

    //antigenic evolution parameters;

    static double antigenicMu_a = 0.0001; //antigenic mutation rate (rate of beneficial mutations)
    static double antigenicMu_b = 0.000005; // antigenic mutation rate for large effect mutations
    static double meanEffect = 0.00038;//0.017; //average selection coefficient (from Sanjuan et al (2004) empirical estimate of mean beneficial mutation mutation effect)
    static double sdEffect = 0.012; //standard deviation of mutational effect;
    //static Gamma dfe = new Gamma((meanEffect*meanEffect)/(sdEffect*sdEffect), meanEffect/(sdEffect*sdEffect), randomGenerator);

    static double s_b1 = meanEffect;
    static double s_b2 = meanEffect;
    static double s_large = 0.014;
    static double s_n1 = 0.005;
    static double s_n2 = 0.005;
    static double e_b = 0.0;
    static double deleteriousMu = 0.000;
    static double neutralMu = 0.0001;
    static double negSel = 0.0001;

    //new 2016-07

    static double rawMutationRate = 0.052;
    static double p_b_seg1 = 0.4;
    static double p_d_seg1 = 0.4;
    static double p_b_seg2 = 0.05;
    static double p_d_seg2 = 0.45;
    static double relativeRate = 50; //relative coinfection rate to infection rate



    static double deleteriousEffect = 0.00000;
    static double mean_n = 0.6;
    static double variance_n = 0.4*0.4;
    static double alpha = mean_n*mean_n/(variance_n);
    static double lambda = mean_n/variance_n;
    static double p = 1; //prob of being in the continuous selection regime
    static double q = 1.0-p; //prob of being in the episodic selection regime

    static double antigenicMu = 0.00075;
    static double nonAntigenicMu = 0.1;
    static double s_ben = 0.2;
    static Exponential benDFE = new Exponential(1/s_ben, randomGenerator);
    static double s_del = 0.1;
    static Exponential delDFE = new Exponential(1/s_del, randomGenerator);
    static double p_ben1 = 0.75;
    static double p_ben2 = 0.15;

    static Normal epistasis = new Normal(0,1, randomGenerator);

    public void print() {

        System.out.println();
        System.out.println("Parameters:");
        System.out.println("...............................................");
 //       System.out.println("pop size, N :"+N);
 //       System.out.println("R0 :"+R0);
 //       System.out.println("Life span (years) :"+lifeSpan);
        System.out.println("Duration of Infection (per day) :"+durationOfInfection);
 //       System.out.println("Waning Immunity (per year) :"+waningImmunity);

        System.out.println("psi (degree of susceptibility reduction) :"+psi);
//        System.out.println("psi_i :"+psi_i);
//        System.out.println("psi_j :"+psi_j);
//        System.out.println("epsilon (strength of seasonality) :"+epsilon);
        System.out.println("Simulation Time (years) :"+simulationTime);

//        System.out.println("Viral Load (constant) :"+viralLoad);
//        System.out.println("Within-host reassortment rate (per day) :"+mutnRate);

//        System.out.println("gamma (per day) :"+gamma);
//        System.out.println("nu_s (per day) :"+nu_s);
//        System.out.println("nu_co (per day) :"+nu_co);
//        System.out.println("mu (per day) :"+mu);
//        System.out.println("beta_s (per day) :"+beta_s);
//        System.out.println("beta_co (per day) :"+beta_co);
        System.out.println();
//        System.out.println("*** Two-patch model parameters ***");
//        System.out.println("epsilon endemic :"+epsilon_endemic);
//        System.out.println("epsilon seasonal :"+epsilon_seasonal);
//        System.out.println("coupling rate, mij (per day) :"+m_ij);

        //initial conditions - start at equilibrium values:
        System.out.println();
        System.out.println("*** Initial conditions ***");
//        System.out.println("S_init :"+S_init);
        System.out.println("Is_init :"+Is_init);
        System.out.println("Ico_init :"+Ico_init);
//        System.out.println("R_init :"+R_init);

//        System.out.println("Is_init_i :"+Is_init_i);
//        System.out.println("Is_init_j :"+Is_init_j);
//        System.out.println("S_init_i :"+S_init_i);
//        System.out.println("S_init_j :"+S_init_j);
//        System.out.println("Ico_init_i :"+Ico_init_i);
//        System.out.println("Ico_init_j :"+Ico_init_j);
//        System.out.println("R_init_i :"+R_init_i);
//        System.out.println("R_init_j :"+R_init_j);

        System.out.println();
        System.out.println("*** Antigenic evolution parameters ***");
//        System.out.println("Antigenic mutation rate1: "+ antigenicMu_a);
//        System.out.println("Antigenic mutation rate2: "+antigenicMu_b);
//        System.out.println("Antigenic mutation rate3: "+ negSel);
//
//        //System.out.println("SD mutational effect: "+sdEffect);
//        System.out.println("Mean dfe: "+s_b1);
//        System.out.println("Large effect: "+s_b2);
//        System.out.println("Del effect seg1: "+s_n1);
//        System.out.println("Del effect seg2: "+s_n2);
//
//
//        System.out.println("Seg2 mutation rate: "+neutralMu);
//        System.out.println("mean: "+mean_n);
//        System.out.println("variance: "+variance_n);



//        System.out.println("p: "+p);
//        System.out.println("q: "+q);

        System.out.println("Antigenic mutation rate: "+antigenicMu);
//        System.out.println("Non-antigenic mutation rate: "+ nonAntigenicMu);
//        System.out.println("s_b: "+s_ben+", s_d: "+s_del);
//        System.out.println("p_b1: "+p_ben1+", p_b2: "+p_ben2);
////        System.out.println("deleterious mutation rate: "+ deleteriousMu);
//        System.out.println("deleterious mutation s: "+ deleteriousEffect);

        System.out.println("p_b_s1: "+p_b_seg1);
        System.out.println("p_b_s2: "+p_b_seg2);
        System.out.println("p_d_s1: "+p_d_seg1);
        System.out.println("p_d_s2: "+p_d_seg2);
        System.out.println("s_ben: "+s_ben);
        System.out.println("s_del: "+s_del);
        System.out.println("mutn_rate: "+rawMutationRate);

        System.out.println("...............................................");


    }

}
