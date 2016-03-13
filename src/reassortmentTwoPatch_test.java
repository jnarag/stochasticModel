/**
 * Created with IntelliJ IDEA.
 * User: Jayna
 * Date: 18/10/2012
 * Time: 14:21
 * To change this template use File | Settings | File Templates.
 */
import cern.colt.list.DoubleArrayList;
import cern.jet.random.Poisson;
import cern.jet.random.engine.RandomEngine;
import cern.jet.random.engine.*;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;
import java.util.Random;

/**
 * Created with IntelliJ IDEA.
 * User: jr207
 * Date: 9/7/12
 * Time: 12:32 PM
 * To change this template use File | Settings | File Templates.
 */
public class reassortmentTwoPatch_test implements EpiModel{


    private double tau = 5.0/16.0;
    public List<int[]> I_matrix = new ArrayList<int[]>();
    public List<int[]> Ico_matrix = new ArrayList<int[]>();
    private RandomEngine randomGenerator;


    @Override
    public void runSimulation(EpiParams params) {

        double maxTime = params.simulationTime*365.25;

        int N_i = params.N;
        //initial numbers;
        int X_curr_i = params.S_init_i;
        int Y_curr_i = params.Is_init_i;
        int Yco_curr_i = params.Ico_init_i;
        int Z_curr_i = params.R_init_i;

        int N_j = params.N;
        //initial numbers;
        int X_curr_j = params.S_init_j;
        int Y_curr_j = params.Is_init_j;
        int Yco_curr_j = params.Ico_init_j;
        int Z_curr_j = params.R_init_j;

        DoubleArrayList rates = new DoubleArrayList();
        double t_curr = 0;
        double t_next;


        double beta_s_i;
        double beta_s_j;
        double beta_co_i;
        double beta_co_j;

        // rates

        double muN_i;     //birth
        double BSIs_i;    //transmission from Is to S
        double BSIco_i;   //transmission from Ico to S
        double muS_i;     //death of S
        double BIsIs_i;   //transmission from Is to Is - coinfection rate
        double nuIs_i;    //recovery of Is
        double nuIco_i;   //recovery of Ico (recovers from both strains simultaneously)
        double muIs_i;    //death of Is
        double muIco_i;   //death of Ico
        double muR_i;     //death of R
        double gammaR_i;  //waning immunity rate, rate at which R flows into S

        double muN_j;     //birth
        double BSIs_j;    //transmission from Is to S
        double BSIco_j;   //transmission from Ico to S
        double muS_j;     //death of S
        double BIsIs_j;   //transmission from Is to Is - coinfection rate
        double nuIs_j;    //recovery of Is
        double nuIco_j;   //recovery of Ico (recovers from both strains simultaneously)
        double muIs_j;    //death of Is
        double muIco_j;   //death of Ico
        double muR_j;     //death of R
        double gammaR_j;  //waning immunity rate, rate at which R flows into S

        double mij_Is_i;
        double mij_Ico_i;
        double mij_S_i;
        double mij_R_i;
        double mij_Is_j;
        double mij_Ico_j;
        double mij_S_j;
        double mij_R_j;


        File outputfile = new File("two_patch_model_reassortants_m_0.002_TEST_psi_i_0.0_psi_j_0.9.csv");
        File iMatrix_output = new File("I_matrix_"+params.m_ij+"_psii_"+params.psi_i+"_"+params.psi_j+"_simTime_"+params.simulationEndTime+".csv");
        File icoMatrix_outout = new File("Ico_matrix_"+params.m_ij+"_psii_"+params.psi_i+"_"+params.psi_j+"_simTime_"+params.simulationEndTime+".csv");

        FileWriter writer = null;
        FileWriter writer1 = null;
        FileWriter writer2 = null;
        try {
            writer = new FileWriter(outputfile);
        } catch (IOException e) {
            e.printStackTrace();  //To change body of catch statement use File | Settings | File Templates.
        }
        try {
            writer1 = new FileWriter(outputfile);
        } catch (IOException e) {
            e.printStackTrace();  //To change body of catch statement use File | Settings | File Templates.
        }
        try {
            writer2 = new FileWriter(outputfile);
        } catch (IOException e) {
            e.printStackTrace();  //To change body of catch statement use File | Settings | File Templates.
        }


        while (t_curr <= maxTime) {


            try {
                writer.write(t_curr+","+X_curr_i+","+Y_curr_i+","+Yco_curr_i+","+Z_curr_i+","+X_curr_j+","+Y_curr_j+","+Yco_curr_j+","+Z_curr_j+"\n");
                writer.flush();
            } catch (IOException e) {
                e.printStackTrace();  //To change body of catch statement use File | Settings | File Templates.
            }

            System.out.println("1 "+X_curr_i+":"+Y_curr_i+":"+Yco_curr_i+":"+Z_curr_i);
            System.out.println("2 "+X_curr_j+":"+Y_curr_j+":"+Yco_curr_j+":"+Z_curr_j);
            rates.removeAll(rates);

            beta_s_i =  params.beta_s*(1.0+params.epsilon_endemic*Math.sin((2.0*Math.PI/365.0)*t_curr));
            beta_s_j = params.beta_s*(1.0+params.epsilon_seasonal*Math.sin((2.0*Math.PI/365.0)*t_curr));

            beta_co_i = beta_s_i;
            beta_co_j = beta_s_j;

            muN_i = params.mu*N_i;
            BSIs_i =beta_s_i*X_curr_i*Y_curr_i/N_i;
            BSIco_i = beta_s_i*X_curr_i*Yco_curr_i/N_i;
            muS_i = params.mu*X_curr_i;
            BIsIs_i = params.psi_i*beta_co_i*Y_curr_i*Y_curr_i/N_i;
            nuIs_i = params.nu_s*Y_curr_i;
            nuIco_i = params.nu_co*Yco_curr_i;
            muIs_i = params.mu*Y_curr_i;
            muIco_i = params.mu*Yco_curr_i;
            muR_i = params.mu*Z_curr_i;
            gammaR_i = params.gamma*Z_curr_i;

            muN_j = params.mu*N_j;
            BSIs_j = beta_s_j*X_curr_j*Y_curr_j/N_j;
            BSIco_j = beta_s_j*X_curr_j*Yco_curr_j/N_j;
            muS_j = params.mu*X_curr_j;
            BIsIs_j = params.psi_j*beta_co_j*Y_curr_j*Y_curr_j/N_j;
            nuIs_j = params.nu_s*Y_curr_j;
            nuIco_j = params.nu_co*Yco_curr_j;
            muIs_j = params.mu*Y_curr_j;
            muIco_j = params.mu*Yco_curr_j;
            muR_j = params.mu*Z_curr_j;
            gammaR_j = params.gamma*Z_curr_j;

            mij_Ico_i = params.m_ij*Yco_curr_i;
            mij_Is_i = params.m_ij*Y_curr_i;
            mij_S_i = params.m_ij*X_curr_i;
            mij_R_i = params.m_ij*Z_curr_i;

            mij_Ico_j = params.m_ij*Yco_curr_j;
            mij_Is_j = params.m_ij*Y_curr_j;
            mij_S_j = params.m_ij*X_curr_j;
            mij_R_j = params.m_ij*Z_curr_j;

            rates.add(muN_i);
            rates.add(BSIs_i);
            rates.add(BSIco_i);
            rates.add(muS_i);
            rates.add(BIsIs_i);
            rates.add(nuIs_i);
            rates.add(nuIco_i);
            rates.add(muIs_i);
            rates.add(muIco_i);
            rates.add(muR_i);
            rates.add(gammaR_i);

            rates.add(muN_j);
            rates.add(BSIs_j);
            rates.add(BSIco_j);
            rates.add(muS_j);
            rates.add(BIsIs_j);
            rates.add(nuIs_j);
            rates.add(nuIco_j);
            rates.add(muIs_j);
            rates.add(muIco_j);
            rates.add(muR_j);
            rates.add(gammaR_j);

            rates.add(mij_Ico_i);
            rates.add(mij_Is_i);
            rates.add(mij_S_i);
            rates.add(mij_R_i);
            rates.add(mij_Ico_j);
            rates.add(mij_Is_j);
            rates.add(mij_S_j);
            rates.add(mij_R_j);

            t_next = tau+t_curr;
            Poisson poisson;


            int event=0;
            int noOfRates = rates.size();
            int minNo;
            int num;
            Random randomGen = new Random();
            while(event<noOfRates) {

                randomGenerator = new MersenneTwister(randomGen.nextInt());
                poisson = new Poisson(tau*rates.get(event), randomGenerator);
                num = poisson.nextInt();
                switch(event) {

                    case 0:  // birth in patch i

                        X_curr_i += num;

                        break;
                    case 1: // transmission in patch i from Is

                        if(Y_curr_i==0){
                            break;
                        }
                        minNo = min(X_curr_i, num);

                        Y_curr_i += minNo;
                        X_curr_i -= minNo;


                        break;
                    case 2:       // transmission from Ico in patch i

                        if(Yco_curr_i==0){
                            break;
                        }
                        minNo = min(X_curr_i, num);

                        Y_curr_i += minNo;
                        X_curr_i -= minNo;

                        break;
                    case 3:       // death of a susceptible in patch i

                        minNo = min(X_curr_i, num);

                        X_curr_i -= minNo;

                        break;
                    case 4:      // coinfection in patch i

                        minNo = min(Y_curr_i, num);

                        Yco_curr_i+= minNo;
                        Y_curr_i-= minNo;

                    case 5:     // recovery in patch i

                        minNo = min(Y_curr_i, num);
                        Y_curr_i-= minNo;
                        Z_curr_i+= minNo;

                        break;
                    case 6:     //recovery in patch i (Ico)

                        minNo = min(Yco_curr_i, num);

                        Yco_curr_i -= minNo;
                        Z_curr_i += minNo;

                        break;
                    case 7:    // death of Is in patch i

                        minNo = min(Y_curr_i, num);
                        Y_curr_i -= minNo;

                        break;
                    case 8:   // death of Ico in patch i

                        minNo = min(Yco_curr_i, num);

                        Yco_curr_i -= minNo;

                        break;
                    case 9:  // death of R in patch i

                        minNo = min(Z_curr_i, num);
                        Z_curr_i -= minNo;

                        break;
                    case 10:  // waning immunity of R in patch i

                        minNo = min(Z_curr_i, num);

                        Z_curr_i -= minNo;
                        X_curr_i += minNo;

                        break;
                    case 11:  // birth in patch j

                        X_curr_j += num;

                        break;
                    case 12: // transmission in patch j from Is

                        if(Y_curr_j==0){
                            break;
                        }
                        minNo = min(X_curr_j, num);
                        Y_curr_j += minNo;
                        X_curr_j -= minNo;

                        break;
                    case 13:       // transmission from Ico in patch j

                        if(Yco_curr_j==0){
                            break;
                        }
                        minNo = min(X_curr_j, num);

                        Y_curr_j += minNo;
                        X_curr_j -= minNo;

                        break;
                    case 14:       // death of a susceptible in patch j

                        minNo = min(X_curr_j, num);

                        X_curr_j -= minNo;

                        break;
                    case 15:      // coinfection in patch j

                        minNo = min(Y_curr_j, num);

                        Yco_curr_j += minNo;
                        Y_curr_j -= minNo;

                        break;
                    case 16:    // recovery of Is in patch j

                        minNo = min(Y_curr_j, num);

                        Y_curr_j -= minNo;
                        Z_curr_j += minNo;
                        break;
                    case 17:   // recovery of Ico in patch j

                        minNo = min(Yco_curr_j, num);

                        Yco_curr_j -= minNo;
                        Z_curr_j += minNo;
                        break;
                    case 18:  // death of Is in patch j

                        minNo = min(Y_curr_j, num);

                        Y_curr_j -= minNo;

                        break;
                    case 19:   // death of Ico in patch j

                        minNo = min(Yco_curr_j, num);

                        Yco_curr_j -= minNo;

                        break;
                    case 20:  // death of R in patch j

                        minNo = min(Z_curr_j, num);

                        Z_curr_j -= minNo;

                        break;
                    case 21:   // waning immunity of R in j

                        minNo = min(Z_curr_j, num);

                        Z_curr_j -= minNo;
                        X_curr_j += minNo;

                        break;
                    case 22:  // migration of Ico from i to j

                        minNo = min(Yco_curr_i, num);

                        Yco_curr_j += minNo;
                        Yco_curr_i -= minNo;

                        break;
                    case 23: // migration of Is from i to j

                        minNo = min(Y_curr_i, num);

                        Y_curr_j += minNo;
                        Y_curr_i -= minNo;

                        break;
                    case 24: // migration of S from i to j

                        minNo = min(X_curr_i, num);

                        X_curr_j += minNo;
                        X_curr_i -= minNo;

                        break;
                    case 25:  // migration of R from i to j

                        minNo = min(Z_curr_i, num);

                        Z_curr_j += minNo;
                        Z_curr_i -= minNo;

                        break;
                    case 26:  // migration of Ico from j to i

                        minNo = min(Yco_curr_j, num);

                        Yco_curr_i += minNo;
                        Yco_curr_j -= minNo;
                        break;
                    case 27: // migration of Is from j to i

                        minNo = min(Y_curr_j, num);

                        Y_curr_i += minNo;
                        Y_curr_j -= minNo;

                        break;
                    case 28:  // migration of S from j to i

                        minNo = min(X_curr_j, num);

                        X_curr_i += minNo;
                        X_curr_j -= minNo;

                        break;
                    case 29:  // migration of R from j to i

                        minNo = min(Z_curr_j, num);

                        Z_curr_i += minNo;
                        Z_curr_j -= minNo;

                        break;

                }
                event++;
            }
            t_curr = t_next;
        }
        try {
            writer.close();
        } catch (IOException e) {
            e.printStackTrace();  //To change body of catch statement use File | Settings | File Templates.
        }


    }

    @Override
    public double getRatesSum(DoubleArrayList rates) {
        double sum = 0;
        for(int i = 0; i < rates.size(); i++) {
            sum += rates.get(i);
        }
        return sum;
    }

    @Override
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

    @Override
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




    public void initializeI_matrix(List<int[]> I_matrix, int I_curr) {
        for(int i=0;i<I_curr;i++) {
            int[] x = new int[6];
            x[0] = -1;
            x[1] = -1;
            x[2] = -1;
            x[3] = -1;
            x[4] = 0;
            x[5] = -1;

            I_matrix.add(x);

        }

    }

    public void initializeIco_matrix(List<int[]> Ico_matrix, int Ico_curr, int I_curr) {
        for(int i=0;i<Ico_curr;i++) {
            int[] x = new int[5];
            x[0] = -1;
            x[1] = -1;
            x[2] = (int)Math.floor(Math.random()*I_curr);
            x[3] = (int)Math.floor(Math.random()*I_curr);
            x[4] = -1;

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
    private int min(int a, int b) {


        if(a>b) {

            return b;
        }
        else {
            return a;
        }
    }
    public infectionHistory getIMatrix() {
        return iMatrix;
    }

    public coinfectionHistory getIcoMatrix() {
        return icoMatrix;
    }
    public static void main(String args[]) {

        EpiParams params = new EpiParams();
        reassortmentTwoPatch_test twoPatchModel = new reassortmentTwoPatch_test();
        twoPatchModel.runSimulation(params);




    }
}

