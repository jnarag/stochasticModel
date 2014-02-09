import cern.jet.random.Exponential;
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
    static int N = (int)7e6;
    static double R0 = 3;
    static int lifeSpan = 30;  //in years
    static int durationOfInfection = 5; // in days
    static int waningImmunity = 4; //in years

    static double psi = 0.0;
    static double psi_i = 1.0;
    static double psi_j = 1.0;
    static double epsilon = 0.0;
    static double simulationTime = 40;
    static double simulationStartTime = 0.0;
    static double simulationEndTime = simulationTime*365.25;

    static int viralLoad = 1000;
    static int mutnRate = 500;

    static double gamma = 1.0/((double)waningImmunity*365.25);
    static double mig = 0.0;
    static double nu_s = 1.0/((double)durationOfInfection);
    static double nu_co = nu_s;
    static double mu = 1.0/((double)lifeSpan*365.25);
    static double beta_s = R0/((double)durationOfInfection);
    static double beta_co = beta_s;

    //initial conditions - start at equilibrium values:
    static int S_init = (int)Math.floor(N/R0);
    static int Is_init = 25000;//(int)Math.floor((mu*N)/(beta_s*(R0-1.0)));
    static int Ico_init = 40;//(int)Math.floor(beta_co/((N*(nu_co+mu))*(Math.pow(Is_init,2))));
    static int R_init = N-S_init-Is_init-Ico_init;

    static double epsilon_endemic = 0.0;
    static double epsilon_seasonal = 0.15;
    static double m_ij = 0.001;
    static int Is_init_i = 25000;
    static int Is_init_j = 25000;
    static int S_init_i = S_init;
    static int S_init_j = S_init;
    static int Ico_init_i = 40;
    static int Ico_init_j = 40;
    static int R_init_i = N-S_init_i-Is_init_i-Ico_init_i;
    static int R_init_j = N-S_init_j-Is_init_j-Ico_init_j;

    static RandomEngine randomGenerator = new MersenneTwister();

    //antigenic evolution parameters;

    static double antigenMu = 0.00; //antigenic mutation rate (rate of beneficial mutations)
    static double meanEffect = 0.001;//0.017; //average selection coefficient (from Sanjuan et al (2004) empirical estimate of mean beneficial mutation mutation effect)
    static double sdEffect = 0.012; //standard deviation of mutational effect;
    //static Gamma dfe = new Gamma((meanEffect*meanEffect)/(sdEffect*sdEffect), meanEffect/(sdEffect*sdEffect), randomGenerator);

    static double dfe = meanEffect;
    static Exponential expDFE = new Exponential(dfe, randomGenerator);


    public void print() {

        System.out.println();
        System.out.println("Parameters:");
        System.out.println("...............................................");
        System.out.println("pop size, N :"+N);
        System.out.println("R0 :"+R0);
        System.out.println("Life span (years) :"+lifeSpan);
        System.out.println("Duration of Infection (per day) :"+durationOfInfection);
        System.out.println("Waning Immunity (per year) :"+waningImmunity);

        System.out.println("psi (degree of susceptibility reduction) :"+psi);
        System.out.println("psi_i :"+psi_i);
        System.out.println("psi_j :"+psi_j);
        System.out.println("epsilon (strength of seasonality) :"+epsilon);
        System.out.println("Simulation Time (years) :"+simulationTime);

        System.out.println("Viral Load (constant) :"+viralLoad);
        System.out.println("Within-host reassortment rate (per day) :"+mutnRate);

        System.out.println("gamma (per day) :"+gamma);
        System.out.println("nu_s (per day) :"+nu_s);
        System.out.println("nu_co (per day) :"+nu_co);
        System.out.println("mu (per day) :"+mu);
        System.out.println("beta_s (per day) :"+beta_s);
        System.out.println("beta_co (per day) :"+beta_co);
        System.out.println();
        System.out.println("*** Two-patch model parameters ***");
        System.out.println("epsilon endemic :"+epsilon_endemic);
        System.out.println("epsilon seasonal :"+epsilon_seasonal);
        System.out.println("coupling rate, mij (per day) :"+m_ij);

        //initial conditions - start at equilibrium values:
        System.out.println();
        System.out.println("*** Initial conditions ***");
        System.out.println("S_init :"+S_init);
        System.out.println("Is_init :"+Is_init);
        System.out.println("Ico_init :"+Ico_init);
        System.out.println("R_init :"+R_init);

        System.out.println("Is_init_i :"+Is_init_i);
        System.out.println("Is_init_j :"+Is_init_j);
        System.out.println("S_init_i :"+S_init_i);
        System.out.println("S_init_j :"+S_init_j);
        System.out.println("Ico_init_i :"+Ico_init_i);
        System.out.println("Ico_init_j :"+Ico_init_j);
        System.out.println("R_init_i :"+R_init_i);
        System.out.println("R_init_j :"+R_init_j);

        System.out.println();
        System.out.println("*** Antigenic evolution parameters ***");
        System.out.println("Antigenic mutation rate: "+ antigenMu);
        //System.out.println("Mean mutational effect: "+meanEffect);
        //System.out.println("SD mutational effect: "+sdEffect);
        System.out.println("Mean dfe: "+dfe);

        System.out.println("...............................................");


    }

}
