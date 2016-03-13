import cern.colt.list.DoubleArrayList;
import cern.jet.random.Exponential;
import cern.jet.random.engine.MersenneTwister;
import cern.jet.random.engine.RandomEngine;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

/**
 * Created with IntelliJ IDEA.
 * User: jr207
 * Date: 8/28/12
 * Time: 11:04 AM
 *
 * Stochastic version of the coinfection influenza model using Gillespie Stochastic Algorithm
 */
public class reassortmentGSA implements EpiModel{

    private RandomEngine randomGenerator = new MersenneTwister();
    public List<int[]> I_matrix = new ArrayList<int[]>();
    public List<int[]> Ico_matrix = new ArrayList<int[]>();



    public void runSimulation(EpiParams params) {


        double maxTime = params.simulationTime*365.25;

        int N = params.N;
        //initial numbers;
        int X_curr = params.S_init;
        int Y_curr = params.Is_init;
        int Yco_curr = params.Ico_init;
        int Z_curr = params.R_init;


        DoubleArrayList rates = new DoubleArrayList();
        double t_curr = 0;
        double t_next;

        // rates

        double muN;     //birth
        double BSIs;    //transmission from Is to S
        double BSIco;   //transmission from Ico to S
        double muS;     //death of S
        double BIsIs;   //transmission from Is to Is - coinfection rate
        double nuIs;    //recovery of Is
        double nuIco;   //recovery of Ico (recovers from both strains simultaneously)
        double muIs;    //death of Is
        double muIco;   //death of Ico
        double muR;     //death of R
        double gammaR;  //waning immunity rate, rate at which R flows into S
        double antigenMu;



        List<Integer> I_matrix_curr = getCurrentInfected(Y_curr, 0);
        List<Integer> Ico_matrix_curr = getCurrentInfected(Yco_curr, 0);
        Map<Integer,Integer> coinfectList = new HashMap<Integer,Integer>();
        int I_matrix_length = (int)Y_curr;
        int Ico_matrix_length = (int)Yco_curr;


        initializeI_matrix(I_matrix,(int)Y_curr);
        initializeIco_matrix(Ico_matrix, (int) Yco_curr, (int) Y_curr);

        //parameters to draw the timestep and choose which event occurs
        Exponential expDist;
        double r1;
        double deltaT;
        double random;
        double rates_sum;
        int event;

        //parameters for the within-host model for coinfected transmission

        double p1;
        double p2;
        int virus;
        double rho;
        double tau;
        double nTrials;
        double reassortantLost;
        double reassortantNotFixedInTau;



        while (t_curr <= maxTime)    {

            System.out.println(t_curr+": " +X_curr+":"+Y_curr+":"+Yco_curr+":"+Z_curr);

            //System.out.println(I_matrix_curr.size()+":"+Ico_matrix_curr.size());
            rates.removeAll(rates);

            muN = params.mu*N;
            BSIs = params.beta_s*X_curr*Y_curr/N;
            BSIco = params.beta_s*X_curr*Yco_curr/N;
            muS = params.mu*X_curr;
            BIsIs = params.psi*params.beta_co*Y_curr*Y_curr/N;
            nuIs = params.nu_s*Y_curr;
            nuIco = params.nu_co*Yco_curr;
            muIs = params.mu*Y_curr;
            muIco = params.mu*Yco_curr;
            muR = params.mu*Z_curr;
            gammaR = params.gamma*Z_curr;

            rates.add(muN);
            rates.add(BSIs);
            rates.add(BSIco);
            rates.add(muS);
            rates.add(BIsIs);
            rates.add(nuIs);
            rates.add(nuIco);
            rates.add(muIs);
            rates.add(muIco);
            rates.add(muR);
            rates.add(gammaR);

            //the events/rates are exponentially distributed with mean = total sum of rates
            rates_sum = getRatesSum(rates);
            expDist = new Exponential(rates_sum, randomGenerator);
            deltaT =  expDist.nextDouble();
            t_next = t_curr + deltaT;
            r1 = randomGenerator.nextDouble();
            int time = (int)Math.floor(t_curr*24*60*60);    //time in sec

            event = chooseEvent(rates, r1, rates_sum);


            //fields in the I_matrix and Ico_matrix
            int[] infectedHistory = new int[5];
            int[] currentInfectedHistory = new int[5];
            int[] coinfectHistory = new int[4];
            int[] currentCoinfectedHistory = new int[4];
            int birth;
            int death;
            int parent_1;
            int parent_2;
            int parent_co;
            int reassortant;

            switch(event) {

                case 0:     // birth of S

                    X_curr++;

                    break;
                case 1:     // transmission from Is

                    Y_curr++;
                    X_curr--;

                    int I_index = (int)Math.floor(Math.random()*I_matrix_curr.size());
                    Integer parent = I_matrix_curr.get(I_index);

                    getInfectedHistory(infectedHistory, time, -1, parent, -1, 0);
                    I_matrix.add(infectedHistory);

                    I_matrix_curr.add(I_matrix_length);
                    I_matrix_length++;

                    break;
                case 2:     // transmission from Ico

                    int Ico_index =  (int)Math.floor(Math.random()*Ico_matrix_curr.size());//currentlyInfected_indices.nextInt();
                    int coinfected = Ico_matrix_curr.get(Ico_index);

                    currentCoinfectedHistory = Ico_matrix.get(coinfected);


                    birth = time;
                    death = -1;
                    parent_1 = -1;     //we can figure this out after we have finished the simulation
                    parent_co = coinfected;


                    //how long has coinfected host been transmitting? Thus, tau = birth-t_curr and we're not using an average tau i.e 1/(Bco*S/N)

                    tau = (currentCoinfectedHistory[0]-time)/(24*60*60);  //So it is scaled in days
                    nTrials = params.mutnRate*tau;
                    reassortantLost = 1.0-1.0/(double)params.viralLoad;
                    reassortantNotFixedInTau = Math.pow(reassortantLost,nTrials);
                    rho = 1.0 - reassortantNotFixedInTau;

                    DoubleArrayList coinfectedTransmission = new DoubleArrayList();

                    coinfectedTransmission.add(rho);
                    coinfectedTransmission.add(reassortantNotFixedInTau/2);
                    coinfectedTransmission.add(reassortantNotFixedInTau/2);

                    random = Math.random();
                    virus = chooseEvent(coinfectedTransmission, random, 1.0);

                    if(virus>0) {

                        reassortant = 0;
                        parent_1 = virus;

                    }
                    else{

                        reassortant = 1;
                    }

                    getInfectedHistory(infectedHistory, birth, death, parent_1, parent_co, reassortant);

                    I_matrix.add(infectedHistory);
                    //new infected

                    Y_curr++;
                    X_curr--;

                    I_matrix_curr.add(I_matrix_length);
                    I_matrix_length++;


                    //for the I_matrix need to see if it is a reassortant or not, and check the coinfected matrix;

                    break;
                case 3:     // death of S

                    X_curr--;

                    break;
                case 4:     // Is contacting another Is transmission - generating a coinfected individual

                    //donor virus/host
                    int p1_index = (int)Math.floor(Math.random()*I_matrix_curr.size());
                    parent_1 = I_matrix_curr.get(p1_index);
                    //I_matrix_curr.remove(p1_index);

                    //recipient virus/host
                    int p2_index = (int)Math.floor(Math.random()*I_matrix_curr.size());
                    parent_2 = I_matrix_curr.get(p2_index);
                    I_matrix_curr.remove(p2_index);

                    birth = time;
                    death = -1;

                    getCoinfectedHistory(coinfectHistory, birth, death, parent_1, parent_2);

                    Ico_matrix.add(coinfectHistory);
                    Ico_matrix_curr.add(Ico_matrix_length);

                    Yco_curr++;
                    Y_curr--;
                    Ico_matrix_length++;

                    break;
                case 5:     // recovery of Is

                    int nuIs_index = (int)Math.floor(Math.random()*I_matrix_curr.size());
                    Integer recovered_Is = I_matrix_curr.get(nuIs_index);
                    I_matrix_curr.remove(nuIs_index);

                    currentInfectedHistory = I_matrix.get(recovered_Is);

                    getInfectedHistory(infectedHistory, currentInfectedHistory[0], time, currentInfectedHistory[2], currentInfectedHistory[3], currentInfectedHistory[4]);

                    I_matrix.set(recovered_Is, infectedHistory);

                    Y_curr--;
                    Z_curr++;

                    break;
                case 6:     // recovery of Ico

                    int nuIco_index = (int)Math.floor(Math.random()*Ico_matrix_curr.size());
                    Integer recovered_Ico = Ico_matrix_curr.get(nuIco_index);
                    Ico_matrix_curr.remove(nuIco_index);

                    currentCoinfectedHistory = Ico_matrix.get(recovered_Ico);
                    getCoinfectedHistory(coinfectHistory, currentCoinfectedHistory[0], time, currentCoinfectedHistory[2], currentCoinfectedHistory[3]);

                    Ico_matrix.set(recovered_Ico, coinfectHistory);

                    Yco_curr--;
                    Z_curr++;

                    break;
                case 7:     // death of Is

                    int muIs_index = (int)Math.floor(Math.random()*I_matrix_curr.size());
                    Integer dead_Is = I_matrix_curr.get(muIs_index);
                    I_matrix_curr.remove(muIs_index);


                    currentInfectedHistory = I_matrix.get(dead_Is);

                    getInfectedHistory(infectedHistory, currentInfectedHistory[0], time, currentInfectedHistory[2], currentInfectedHistory[3], currentInfectedHistory[4]);

                    I_matrix.set(dead_Is, infectedHistory);

                    Y_curr--;

                    break;
                case 8:     // death of Ico

                    int muIco_index = (int)Math.floor(Math.random()*Ico_matrix_curr.size());
                    Integer dead_Ico = Ico_matrix_curr.get(muIco_index);
                    Ico_matrix_curr.remove(muIco_index);

                    coinfectHistory = Ico_matrix.get(dead_Ico);

                    currentCoinfectedHistory = Ico_matrix.get(dead_Ico);
                    getCoinfectedHistory(coinfectHistory, currentCoinfectedHistory[0], time, currentCoinfectedHistory[2], currentCoinfectedHistory[3]);

                    Ico_matrix.set(dead_Ico, coinfectHistory);

                    Yco_curr--;

                    break;
                case 9:     // death of R

                    Z_curr--;

                    break;
                case 10:    // rate of waning immunity of R individual

                    Z_curr--;
                    X_curr++;

                    break;
            }

            t_curr = t_next;

        }


    }

    public double getRatesSum(DoubleArrayList rates) {

        double sum = 0;
        for(int i = 0; i < rates.size(); i++) {
            sum += rates.get(i);
        }
        return sum;
    }

    public int chooseEvent(DoubleArrayList rates, double randomNo, double rates_sum) {

        double[] cumsum = normCumSum(rates, rates_sum);
        int event = 0;

        while(cumsum[event] < randomNo) {

            event++;

            if(event>=cumsum.length){
                break;
            }

        }

        return event;

    }

    public double[] normCumSum(DoubleArrayList rates, double rateMatrix_sum) {

        double[] cumsum = new double[rates.size()];
        int j = 0;
        for(int i=0;i<rates.size(); i++) {

            cumsum[i] = rates.get(i);

            int k = 0;
            while(k<j) {

                cumsum[i] += rates.get(k);
                k++;

            }
            cumsum[i]/=rateMatrix_sum;
            j++;
        }

        return cumsum;
    }

    private void getInfectedHistory(int[] infectedHistory, int birth, int death, int parent, int parent_co, int reassortant) {

        infectedHistory[0] = birth;
        infectedHistory[1] = death;
        infectedHistory[2] = parent;
        infectedHistory[3] = parent_co;
        infectedHistory[4] = reassortant;

    }

    public void initializeI_matrix(List<int[]> I_matrix, int I_curr) {

        for(int i=0;i<I_curr;i++) {
            int[] x = new int[5];
            x[0] = -1;
            x[1] = -1;
            x[2] = -1;
            x[3] = -1;
            x[4] = -1;

            I_matrix.add(x);

        }
    }

    public void initializeIco_matrix(List<int[]> Ico_matrix, int Ico_curr, int I_curr) {
        int[] x;
        for(int i=0;i<Ico_curr;i++) {
            x = new int[4];
            x[0] = -1;
            x[1] = -1;
            x[2] = (int)Math.floor(Math.random()*I_curr);
            x[3] = (int)Math.floor(Math.random()*I_curr);

            Ico_matrix.add(x);
        }
    }

    public List<Integer> getCurrentInfected(int I_curr, int startingNumber) {

        List<Integer> I_matrix_curr = new ArrayList<Integer>();
        int i = startingNumber;
        while(i<I_curr){

            I_matrix_curr.add(i);
            i++;
        }

        return I_matrix_curr;
    }

    public infectionHistory getIMatrix() {
        return iMatrix;
    }

    public coinfectionHistory getIcoMatrix() {
        return icoMatrix;
    }

    private void getCoinfectedHistory(int[] coinfectHistory, int birth, int death, int parent1, int parent2) {

        coinfectHistory[0] = birth;
        coinfectHistory[1] = death;
        coinfectHistory[2] = parent1;
        coinfectHistory[3] = parent2;

    }

    public static void main(String [] args) {

        EpiParams params = new EpiParams();
        reassortmentGSA model = new reassortmentGSA();
        model.runSimulation(params);

    }



}
