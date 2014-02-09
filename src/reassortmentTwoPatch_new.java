/**
 * Created with IntelliJ IDEA.
 * User: jr207
 * Date: 11/20/12
 * Time: 5:05 PM
 * To change this template use File | Settings | File Templates.
 */
import cern.colt.list.DoubleArrayList;
import cern.jet.random.Normal;
import cern.jet.random.Poisson;
import cern.jet.random.engine.MersenneTwister;
import cern.jet.random.engine.RandomEngine;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.*;
/**
 * Created with IntelliJ IDEA.
 * User: jr207
 * Date: 9/7/12
 * Time: 12:32 PM
 * To change this template use File | Settings | File Templates.
 */
public class reassortmentTwoPatch_new implements EpiModel{


    private double tau = 1.0/4.0;
    private Random random = new Random();
    private RandomEngine randomGenerator = new MersenneTwister(random.nextInt());
    public infectionHistory iMatrix = new infectionHistory((int)100e6);
    public coinfectionHistory icoMatrix = new coinfectionHistory();
    private FileWriter writer = null;
    private FileWriter writer1 = null;
    private FileWriter writer2 = null;
    private FileWriter writer3 = null;
    private boolean writeOutput = false;
    private boolean writeIncid = true;

    @Override
    public void runSimulation(EpiParams params) {

        params.print();

        icoMatrix.initialize((int)1e6);

        File outputfile = new File("two_patch_model_reassortants_antigenicMu_"+params.antigenMu+"mij_"+params.m_ij+"_psij_"+params.psi_i+"_"+params.psi_j+"_epsij_"+params.epsilon_endemic+"_"+params.epsilon_seasonal+
                "_D_"+params.durationOfInfection+"_W_"+params.waningImmunity+"_R0_"+params.R0+"_simTime_"+params.simulationTime+"yrs.csv");

//        File iMatrix_output = new File("I_matrix_antigenicMu_"+params.antigenMu+"mij_"+params.m_ij+"_psii_"+params.psi_i+"_"+params.psi_j+"_simTime_"+params.simulationTime+"yrs.csv");
//        File icoMatrix_outout = new File("Ico_matrix_mij_"+params.m_ij+"_psii_"+params.psi_i+"_"+params.psi_j+"_simTime_"+params.simulationTime+"yrs.csv");

        File cumIncid = new File("cumulative_Incidence_antigenicMu_"+params.antigenMu+"mij_"+params.m_ij+"_psij_"+params.psi_i+"_"+params.psi_j+"_epsij_"+params.epsilon_endemic+"_"+params.epsilon_seasonal+
                "_D_"+params.durationOfInfection+"_W_"+params.waningImmunity+"_R0_"+params.R0+"_simTime_"+params.simulationTime+"yrs.csv");
        //boolean writeOutput = false;
        try {
            writer = new FileWriter(outputfile);
        } catch (IOException e) {
            e.printStackTrace();  //To change body of catch statement use File | Settings | File Templates.
        }

//        if(writeOutput) {
//            try {
//                writer1 = new FileWriter(iMatrix_output);
//            } catch (IOException e) {
//                e.printStackTrace();  //To change body of catch statement use File | Settings | File Templates.
//            }
//            try {
//                writer2 = new FileWriter(icoMatrix_outout);
//            } catch (IOException e) {
//                e.printStackTrace();  //To change body of catch statement use File | Settings | File Templates.
//            }
//        }

        if(writeIncid){
            try {
                writer3 = new FileWriter(cumIncid);
            } catch (IOException e) {
                e.printStackTrace();  //To change body of catch statement use File | Settings | File Templates.
            }

        }


        double maxTime = params.simulationTime*365.25;

        int N_i = params.N;
        //initial numbers;
        int X_curr_i = params.S_init_i;
        int Y_curr_i = params.Is_init_i;
        int Yco_curr_i = params.Ico_init_i;
        int Z_curr_i = params.R_init_i;
        int Yr_primary_curr_i = 0;
        int Yr_curr_i=0;

        int N_j = params.N;
        //initial numbers;
        int X_curr_j = params.S_init_j;
        int Y_curr_j = params.Is_init_j;
        int Yco_curr_j = params.Ico_init_j;
        int Z_curr_j = params.R_init_j;
        int Yr_primary_curr_j = 0;
        int Yr_curr_j = 0;

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

        double evolve_i; // antigenic evolution in patch i
        double evolve_j; // antigenic evolution in patch j

        List<Integer> I_matrix_curr_i = getCurrentInfected(Y_curr_i, 0);  //current infections in patch i
        List<Integer> Ico_matrix_curr_i = getCurrentInfected(Yco_curr_i, 0); //current coinfections in patch i

        List<Integer> I_matrix_curr_j = getCurrentInfected((Y_curr_i+Y_curr_j), Y_curr_i); //current infections in patch j
        List<Integer> Ico_matrix_curr_j = getCurrentInfected((Yco_curr_i+Yco_curr_j), Yco_curr_i); //current coinfections in patch j

        int I_matrix_length = Y_curr_i+Y_curr_j;
        int Ico_matrix_length = Yco_curr_i+Yco_curr_j;

//        //Fitness parameters;
//        Normal fitnessDistribution = new Normal(6.0, 2.0, randomGenerator);
//        Normal errorDistribution = new Normal(0.0, 2.0, randomGenerator);

        initializeI_matrix(Y_curr_i, Y_curr_j);
        initializeIco_matrix(Yco_curr_i, Yco_curr_j, (Y_curr_i+Y_curr_j));

        //cumulative incidence
        int cumIs_i = 0;
        int cumIco_i = 0;
        int cumIrp_i = 0;
        int cumIrs_i = 0;

        int cumIs_j = 0;
        int cumIco_j = 0;
        int cumIrp_j = 0;
        int cumIrs_j = 0;

        while (t_curr <= maxTime) {


            System.out.println(">"+t_curr+","+X_curr_i+","+Y_curr_i+","+Yco_curr_i+","+Yr_primary_curr_i+","+Yr_curr_i+","+Z_curr_i+","+X_curr_j+","+Y_curr_j+","+Yco_curr_j+","+Yr_primary_curr_j+","+Yr_curr_j+","+Z_curr_j);
            try {
                writer.write(t_curr+","+X_curr_i+","+Y_curr_i+","+Yco_curr_i+","+Yr_primary_curr_i+","+Yr_curr_i+","+Z_curr_i+","+X_curr_j+","+Y_curr_j+","+Yco_curr_j+","+Yr_primary_curr_j+","+Yr_curr_j+","+Z_curr_j+"\n");
                writer.flush();
            } catch (IOException e) {
                e.printStackTrace();  //To change body of catch statement use File | Settings | File Templates.
            }

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

            evolve_i = params.antigenMu*Y_curr_i;     // antigenic evolution in patch 1
            evolve_j = params.antigenMu*Y_curr_j;     // antigenic evolution in patch 2

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

            rates.add(evolve_i);
            rates.add(evolve_j);

            t_next = tau+t_curr;

            int event=0;
            int noOfRates = rates.size();
            int j;

            double birth;
            double death;
            int parent_1;
            int parent_2;
            int parent_co;
            int reassortant;
            int patch;
            int num;
            //Random randomGen = new Random();
            int minNo;

            List<Double> fitnessWeights;

            if(writeIncid) {

                try {
                    writer3.write(t_curr+","+cumIs_i+","+cumIco_i+","+cumIrp_i+","+cumIrs_i+","+cumIs_j+","+cumIco_j+","+cumIrp_j+","+cumIrs_j+"\n");
                    writer3.flush();
                } catch (IOException e) {
                    e.printStackTrace();  //To change body of catch statement use File | Settings | File Templates.
                }

            }
            while(event<noOfRates) {
                Poisson poisson = new Poisson(tau*rates.get(event), randomGenerator);
                num = poisson.nextInt();

                switch(event) {

                    case 0:  // birth in patch i

                        X_curr_i += num;

                        break;
                    case 1: // transmission in patch i from Is

                        if(Y_curr_i==0){
                            break;
                        }
                        minNo = Math.min(X_curr_i, num);

                        //Math.min(X_curr_i, num)

                        Y_curr_i += minNo;
                        X_curr_i -= minNo;
                        patch = 1;

                        //cumulative incid
                        cumIs_i += minNo;
                        j = 0;

                        fitnessWeights = calculateWeights(I_matrix_curr_i);
                        while(j<minNo) {

                            //I want to choose the parent index according to the fitness weights

                            int bin = chooseBin(fitnessWeights);
                            Integer parent = I_matrix_curr_i.get(bin);

                            //find out if the infected host carries a reassortant strain
                            reassortant = iMatrix.getReassortant(parent);//
                            if(reassortant == 1) {
                                Yr_curr_i++;
                                //cumulative Irs_i (secondary Ir)
                                cumIrs_i++;
                            }

                            Double newFitness = iMatrix.getFitness(parent);

                            //log infection in the corresponding vectors;
                            iMatrix.logId(I_matrix_length);
                            iMatrix.logBirth(t_curr);
                            iMatrix.logDeath(Double.NEGATIVE_INFINITY);
                            iMatrix.logParent(parent);
                            iMatrix.logParentCo((int) Double.NEGATIVE_INFINITY);
                            iMatrix.logReassortant(reassortant);
                            iMatrix.logFitness(newFitness);
                            iMatrix.logSegFitness(iMatrix.getSegFitness(parent));
                            iMatrix.logPatch(patch);

                            I_matrix_curr_i.add(I_matrix_length);
                            I_matrix_length++;

                            j++;

                        }
                        break;
                    case 2:       // transmission from Ico in patch i

                        if(Yco_curr_i==0){
                            break;
                        }

                        minNo = Math.min(X_curr_i, num);

                        Y_curr_i += minNo;
                        X_curr_i -= minNo;

                        //cumuative incidence
                        cumIs_i += minNo;

                        patch = 1;
                        j = 0;

                        while(j<minNo) {

                            int Ico_index = (int)Math.floor(Math.random()*Ico_matrix_curr_i.size());
                            int parentCo = Ico_matrix_curr_i.get(Ico_index);

                            int parent = (int)Double.NEGATIVE_INFINITY;

                            //how long has coinfected host been transmitting? Thus, tau = birth-t_curr and we're not using an average tau i.e 1/(Bco*S/N)

                            //double transmTime = t_curr-(Double)icoMatrixBirth.get(parentCo);  //So it is scaled in days
                            double transmTime = t_curr-icoMatrix.getBirth(parentCo);  //So it is scaled in days
                            double nTrials = params.mutnRate*transmTime;
                            double reassortantLost = 1.0-1.0/(double)params.viralLoad;
                            double reassortantNotFixedInTau = Math.pow(reassortantLost,nTrials);
                            double rho = 1.0 - reassortantNotFixedInTau;

                            DoubleArrayList coinfectedTransmission = new DoubleArrayList();

                            coinfectedTransmission.add(rho);
                            coinfectedTransmission.add(reassortantNotFixedInTau/2);
                            coinfectedTransmission.add(reassortantNotFixedInTau/2);

                            double random = Math.random();
                            int virus = chooseEvent(coinfectedTransmission, random, 1.0);

                            double newFitness = 0.0;double seg1Fitness = 0.0;
                            if(virus>0) {

                                reassortant = 0;
                                if(virus==1){

                                    parent = icoMatrix.getParent1(parentCo);//Integer)icoMatrixParent1.get(parentCo);
                                    newFitness = iMatrix.getFitness(parent);//+errorDistribution.nextDouble();
                                    seg1Fitness = iMatrix.getSegFitness(parent);
                                }
                                else{
//
                                    parent = icoMatrix.getParent2(parentCo);//Integer)icoMatrixParent2.get(parentCo);
                                    newFitness = iMatrix.getFitness(parent);//+errorDistribution.nextDouble();
                                    seg1Fitness = iMatrix.getSegFitness(parent);
                                }

                            }
                            else{

                                reassortant = 1;

                                double parent1Fitness = iMatrix.getFitness(icoMatrix.getParent1(parentCo));
                                double parent2Fitness = iMatrix.getFitness(icoMatrix.getParent2(parentCo));

                                double seg1Fitness_p1 = iMatrix.getSegFitness(icoMatrix.getParent1(parentCo));
                                double seg1Fitness_p2 = iMatrix.getSegFitness(icoMatrix.getParent2(parentCo));

                                double seg2Fitness_p1 = parent1Fitness-seg1Fitness_p1;
                                double seg2Fitness_p2 = parent2Fitness-seg1Fitness_p2;

                                int r = (int)Math.round(Math.random());

                                if(r < 1) {  // seg 1 inherited from parent 1 and seg 2 inherited from parent 2


                                    newFitness = seg1Fitness_p1 + seg2Fitness_p2;
                                    seg1Fitness = seg1Fitness_p1;
                                }
                                else{       // seg 2 inherited from parent 2 and seg 2 inherited from parent 1

                                    newFitness = seg1Fitness_p2 + seg2Fitness_p1;
                                    seg1Fitness = seg1Fitness_p2;

                                }

                                //newFitness = (parent1Fitness + parent2Fitness)*0.5 + errorDistribution.nextDouble();
                                //newFitness = ((parent1Fitness + errorDistribution.nextDouble()) + (parent2Fitness + errorDistribution.nextDouble()))*0.5;
                                // newFitness = (parent1Fitness + parent2Fitness)*0.5;
                                Yr_primary_curr_i++;
                                Yr_curr_i++;

                                //cumulative Irp and Irs in patch i
                                cumIrp_i++;
                                cumIrs_i++;

                            }
                            iMatrix.logId(I_matrix_length);
                            iMatrix.logBirth(t_curr);
                            iMatrix.logDeath(Double.NEGATIVE_INFINITY);
                            iMatrix.logParent(parent);
                            iMatrix.logParentCo(parentCo);
                            iMatrix.logReassortant(reassortant);
                            iMatrix.logFitness(newFitness);
                            iMatrix.logSegFitness(seg1Fitness);
                            iMatrix.logPatch(patch);

                           I_matrix_curr_i.add(I_matrix_length);
                            I_matrix_length++;


                            j++;
                        }

                        break;
                    case 3:       // death of a susceptible in patch i

                        minNo = Math.min(X_curr_i, num);
                        X_curr_i -= minNo;

                        break;
                    case 4:      // coinfection in patch i

                        minNo = Math.min(Y_curr_i, num);

                        Yco_curr_i += minNo;
                        Y_curr_i -= minNo;

                        //cumulative Ico incid
                        cumIco_i += minNo;

                        patch = 1;
                        j = 0;

                        while(j<minNo){

                            int recipient_index = (int)Math.floor(Math.random()*I_matrix_curr_i.size());
                            Integer recipientParent = I_matrix_curr_i.get(recipient_index);
                            I_matrix_curr_i.remove(recipientParent);

                            //donor virus/host

                            int donor_index = (int)Math.floor(Math.random()*I_matrix_curr_i.size());
                            int donorParent = I_matrix_curr_i.get(donor_index);

                            reassortant = iMatrix.getReassortant(recipientParent);
                            //should I have update the I_matrix for the recipient virus - will appear as if it is dead
                            if(reassortant == 1) {

                                parent_co = iMatrix.getParentCo(recipientParent);
                                if(parent_co >= 0) {

                                    if(Yr_primary_curr_i>0){
                                        Yr_primary_curr_i--;
                                    }

                                }

                                Yr_curr_i--;
                            }

                            icoMatrix.logBirth(t_curr);
                            icoMatrix.logDeath(Double.NEGATIVE_INFINITY);
                            icoMatrix.logParent1(donorParent);
                            icoMatrix.logParent2(recipientParent);
                            //icoMatrix.logFitness(fitnessDistribution.nextDouble());
                            icoMatrix.logPatch(patch);


                            Ico_matrix_curr_i.add(Ico_matrix_length);
                            Ico_matrix_length++;

                            j++;

                        }
                        break;
                    case 5:   //recovery

                        minNo = Math.min(Y_curr_i, num);
                        Y_curr_i -= minNo;
                        Z_curr_i += minNo;
                        patch = 1;
                        j = 0;
                        while(j<minNo){

                            int nuIs_index = (int)Math.floor(Math.random()*I_matrix_curr_i.size());
                            Integer recovered_Is = I_matrix_curr_i.get(nuIs_index);
                            I_matrix_curr_i.remove(nuIs_index);

                            reassortant = iMatrix.getReassortant(recovered_Is);
                            if(reassortant==1) {

                                parent_co = iMatrix.getParentCo(recovered_Is);
                                if(parent_co >= 0) {
                                    //check if primary reassortant
                                    if(Yr_primary_curr_i>0) {
                                        Yr_primary_curr_i--;
                                    }

                                }

                                Yr_curr_i--;
                            }

                            iMatrix.setDeath(recovered_Is, t_curr);
                            iMatrix.setPatch(recovered_Is, patch);
                            j++;
                        }

                        break;
                    case 6:   //recovery

                        minNo = Math.min(Yco_curr_i, num);

                        Yco_curr_i -= minNo;
                        Z_curr_i += minNo;

                        patch = 1;

                        j = 0;
                        while(j<minNo) {

                            int nuIco_index = (int)Math.floor(Math.random()*Ico_matrix_curr_i.size());
                            Integer recovered_Ico = Ico_matrix_curr_i.get(nuIco_index);
                            Ico_matrix_curr_i.remove(nuIco_index);

                            icoMatrix.setDeath(recovered_Ico, t_curr);
                            icoMatrix.setPatch(recovered_Ico, patch);

                            j++;
                        }

                        break;
                    case 7:   //death

                        minNo = Math.min(Y_curr_i, num);
                        Y_curr_i -= minNo;
                        j = 0;
                        while(j<minNo) {


                            int muIs_index = (int)Math.floor(Math.random()*I_matrix_curr_i.size());
                            Integer dead_Is = I_matrix_curr_i.get(muIs_index);
                            I_matrix_curr_i.remove(muIs_index);

                            reassortant = iMatrix.getReassortant(dead_Is);
                            if(reassortant==1) {

                                parent_co = iMatrix.getParentCo(dead_Is);
                                if(parent_co>=0) {

                                    if(Yr_primary_curr_i>0) {
                                        Yr_primary_curr_i--;
                                    }

                                }

                                Yr_curr_i--;
                            }

                            iMatrix.setDeath(dead_Is, t_curr);


                            j++;
                        }

                        break;
                    case 8:  //death

                        minNo = Math.min(Yco_curr_i, num);

                        Yco_curr_i -= minNo;

                        j = 0;
                        while(j<minNo) {

                            int muIco_index = (int)Math.floor(Math.random()*Ico_matrix_curr_i.size());
                            Integer dead_Ico = Ico_matrix_curr_i.get(muIco_index);
                            Ico_matrix_curr_i.remove(muIco_index);

                            icoMatrix.setDeath(dead_Ico, t_curr);

                            j++;
                        }
                        break;
                    case 9:  //death

                        minNo = Math.min(Z_curr_i, num);
                        Z_curr_i -= minNo;

                        break;
                    case 10: //waning immunity

                        minNo = Math.min(Z_curr_i, num);

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
                        minNo = Math.min(X_curr_j, num);
                        Y_curr_j += minNo;
                        X_curr_j -= minNo;

                        //cumulative Is incidence in patch j
                        cumIs_j += minNo;

                        patch = 2;

                        j = 0;

                        fitnessWeights = calculateWeights(I_matrix_curr_j);
                        while(j<minNo) {

                            int I_index = chooseBin(fitnessWeights);

                            Integer parent = I_matrix_curr_j.get(I_index);
                            reassortant = iMatrix.getReassortant(parent);

                            if(reassortant==1) {
                                Yr_curr_j++;
                                //cumulative incidence Irs in patch j
                                cumIrs_j++;
                            }

                            Double newFitness = iMatrix.getFitness(parent);// +errorDistribution.nextDouble();

                            iMatrix.logId(I_matrix_length);
                            iMatrix.logBirth(t_curr);
                            iMatrix.logDeath(Double.NEGATIVE_INFINITY);
                            iMatrix.logParent(parent);
                            iMatrix.logParentCo((int) Double.NEGATIVE_INFINITY);
                            iMatrix.logReassortant(reassortant);
                            iMatrix.logFitness(newFitness);
                            iMatrix.logSegFitness(iMatrix.getSegFitness(parent));
                            iMatrix.logPatch(patch);

                         I_matrix_curr_j.add(I_matrix_length);
                            I_matrix_length++;

                            j++;
                        }
                        break;
                    case 13:       // transmission from Ico in patch j

                        if(Yco_curr_j==0){
                            break;
                        }
                        minNo = Math.min(X_curr_j, num);

                        Y_curr_j += minNo;
                        X_curr_j -= minNo;

                        //cumulative incidence of Is in patch j
                        cumIs_j += minNo;

                        patch = 2;
                        j = 0;

                        while(j<minNo) {

                            int Ico_index = (int)Math.floor(Math.random()*Ico_matrix_curr_j.size());
                            int parentCo = Ico_matrix_curr_j.get(Ico_index);

                            //currentCoinfectedHistory = Ico_matrix.get(coinfected);

                            int parent = (int)Double.NEGATIVE_INFINITY;     //we can figure this out after we have finished the simulation
                            //parent_co = coinfected;

                            //how long has coinfected host been transmitting? Thus, tau = birth-t_curr and we're not using an average tau i.e 1/(Bco*S/N)

                            //double transmTime = t_curr-(Double)icoMatrixBirth.get(coinfected);  //So it is scaled in days
                            double transmTime = t_curr-icoMatrix.getBirth(parentCo);
                            double nTrials = params.mutnRate*transmTime;
                            double reassortantLost = 1.0-1.0/(double)params.viralLoad;
                            double reassortantNotFixedInTau = Math.pow(reassortantLost,nTrials);
                            double rho = 1.0 - reassortantNotFixedInTau;

                            DoubleArrayList coinfectedTransmission = new DoubleArrayList();

                            coinfectedTransmission.add(rho);
                            coinfectedTransmission.add(reassortantNotFixedInTau/2);
                            coinfectedTransmission.add(reassortantNotFixedInTau/2);

                            double random = Math.random();
                            int virus = chooseEvent(coinfectedTransmission, random, 1.0);
                            double newFitness = 0.0;double seg1Fitness = 0.0;
                            if(virus>0) {

                                reassortant = 0;
                                if(virus==1){

                                    parent = icoMatrix.getParent1(parentCo);
                                    newFitness = iMatrix.getFitness(parent);// + errorDistribution.nextDouble();
                                    seg1Fitness = iMatrix.getSegFitness(parent);
                                }
                                else{

                                    parent = icoMatrix.getParent2(parentCo);
                                    newFitness = iMatrix.getFitness(parent);// + errorDistribution.nextDouble();
                                    seg1Fitness = iMatrix.getSegFitness(parent);
                                }

                            }
                            else{

                                reassortant = 1;

                                double parent1Fitness = iMatrix.getFitness(icoMatrix.getParent1(parentCo));
                                double parent2Fitness = iMatrix.getFitness(icoMatrix.getParent2(parentCo));

                                double seg1Fitness_p1 = iMatrix.getSegFitness(icoMatrix.getParent1(parentCo));
                                double seg1Fitness_p2 = iMatrix.getSegFitness(icoMatrix.getParent2(parentCo));

                                double seg2Fitness_p1 = parent1Fitness-seg1Fitness_p1;
                                double seg2Fitness_p2 = parent2Fitness-seg1Fitness_p2;

                                int r = (int)Math.round(Math.random());
                                //newFitness = (parent1Fitness + parent2Fitness)*0.5 + errorDistribution.nextDouble();
                                //newFitness = ((parent1Fitness + errorDistribution.nextDouble()) + (parent2Fitness + errorDistribution.nextDouble()))*0.5;

                                if(r < 1) {  // seg 1 inherited from parent 1 and seg 2 inherited from parent 2


                                    newFitness = seg1Fitness_p1 + seg2Fitness_p2;
                                    seg1Fitness = seg1Fitness_p1;
                                }
                                else{       // seg 2 inherited from parent 2 and seg 2 inherited from parent 1

                                    newFitness = seg1Fitness_p2 + seg2Fitness_p1;
                                    seg1Fitness = seg1Fitness_p2;

                                }


                                Yr_primary_curr_j++;
                                Yr_curr_j++;

                                //cumulative incidence of primary Ir in patch j
                                cumIrp_j++;
                                cumIrs_j++;
                            }


                            iMatrix.logId(I_matrix_length);
                            iMatrix.logBirth(t_curr);
                            iMatrix.logDeath(Double.NEGATIVE_INFINITY);
                            iMatrix.logParent(parent);
                            iMatrix.logParentCo(parentCo);
                            iMatrix.logReassortant(reassortant);
                            iMatrix.logFitness(newFitness);
                            iMatrix.logSegFitness(seg1Fitness);
                            iMatrix.logPatch(patch);
                            I_matrix_curr_j.add(I_matrix_length);
                            I_matrix_length++;

                            j++;
                        }
                        break;
                    case 14:       // death of a susceptible in patch j

                        minNo = Math.min(X_curr_j, num);

                        X_curr_j -= minNo;

                        break;
                    case 15:      // coinfection in patch j

                        minNo = Math.min(Y_curr_j, num);

                        Yco_curr_j += minNo;
                        Y_curr_j -= minNo;

                        //cumulative incidence
                        cumIco_j += minNo;

                        patch = 2;
                        j = 0;

                        while(j<minNo){


                            //recipient virus/host
                            int recipient_index = (int)Math.floor(Math.random()*I_matrix_curr_j.size());
                            Integer recipientParent = I_matrix_curr_j.get(recipient_index);
                            I_matrix_curr_j.remove(recipientParent);

                            //donor virus/host
                            int donor_index = (int)Math.floor(Math.random()*I_matrix_curr_j.size());
                            Integer donorParent = I_matrix_curr_j.get(donor_index);

                            reassortant = iMatrix.getReassortant(recipientParent);

                            if(reassortant==1) {

                                parent_co = iMatrix.getParentCo(recipientParent);
                                if(parent_co>=0) {

                                    if(Yr_primary_curr_j>0){
                                        Yr_primary_curr_j--;
                                    }
                                }

                                Yr_curr_j--;
                            }

                            icoMatrix.logBirth(t_curr);
                            icoMatrix.logDeath(Double.NEGATIVE_INFINITY);
                            icoMatrix.logParent1(donorParent);
                            icoMatrix.logParent2(recipientParent);
                            icoMatrix.logPatch(patch);
                            //icoMatrix.logFitness(fitnessDistribution.nextDouble());

                            if(writeOutput) {
                                try {
                                    writer2.write(Ico_matrix_length + "," + icoMatrix.getBirth(Ico_matrix_length) + "," + icoMatrix.getDeath(Ico_matrix_length) + "," + icoMatrix.getParent1(Ico_matrix_length) + "," +
                                            icoMatrix.getParent2(Ico_matrix_length)+ "," + icoMatrix.getPatch(Ico_matrix_length)+ "\n");
                                    writer2.flush();
                                } catch (IOException e) {
                                    e.printStackTrace();  //To change body of catch statement use File | Settings | File Templates.
                                }
                            }

                            Ico_matrix_curr_j.add(Ico_matrix_length);
                            Ico_matrix_length++;

                            j++;
                        }
                        break;
                    case 16:

                        //recovery of Is_j

                        minNo = Math.min(Y_curr_j, num);

                        Y_curr_j -= minNo;
                        Z_curr_j += minNo;
                        patch = 2;

                        j = 0;
                        while(j<minNo){

                            int nuIs_index = (int)Math.floor(Math.random()*I_matrix_curr_j.size());
                            Integer recovered_Is = I_matrix_curr_j.get(nuIs_index);
                            I_matrix_curr_j.remove(nuIs_index);

                            reassortant = iMatrix.getReassortant(recovered_Is);
                            if(reassortant==1) {

                                //parent_co = (Integer)iMatrixParent_co.get(recovered_Is);
                                parent_co = iMatrix.getParentCo(recovered_Is);
                                if(parent_co>=0) {

                                    if(Yr_primary_curr_j>0) {
                                        Yr_primary_curr_j--;
                                    }
                                }

                                Yr_curr_j--;

                            }

                            iMatrix.setDeath(recovered_Is, t_curr);
                            iMatrix.setPatch(recovered_Is, patch);

                            j++;
                        }
                        break;
                    case 17:

                        minNo = Math.min(Yco_curr_j, num);

                        Yco_curr_j -= minNo;
                        Z_curr_j += minNo;
                        patch = 2;

                        j = 0;
                        while(j<minNo) {

                            int nuIco_index = (int)Math.floor(Math.random()*Ico_matrix_curr_j.size());
                            Integer recovered_Ico = Ico_matrix_curr_j.get(nuIco_index);
                            Ico_matrix_curr_j.remove(nuIco_index);

                            icoMatrix.setDeath(recovered_Ico, t_curr);
                            icoMatrix.setPatch(recovered_Ico, patch);

                            j++;
                        }
                        break;
                    case 18:
                        //death Is_j
                        minNo = Math.min(Y_curr_j, num);

                        Y_curr_j -= minNo;

                        j = 0;
                        while(j<minNo) {

                            int muIs_index = (int)Math.floor(Math.random()*I_matrix_curr_j.size());
                            Integer dead_Is = I_matrix_curr_j.get(muIs_index);
                            I_matrix_curr_j.remove(muIs_index);

                            reassortant = iMatrix.getReassortant(dead_Is);
                            if(reassortant==1) {

                                parent_co = iMatrix.getParentCo(dead_Is);
                                if(parent_co>=0) {

                                    if(Yr_primary_curr_j>0) {
                                        Yr_primary_curr_j--;
                                    }
                                }

                                Yr_curr_j--;

                            }

                            iMatrix.setDeath(dead_Is, t_curr);

                            j++;
                        }
                        break;
                    case 19:
                        //death Ico_j

                        minNo = Math.min(Yco_curr_j, num);

                        Yco_curr_j -= minNo;

                        j = 0;
                        while(j<minNo) {

                            int muIco_index = (int)Math.floor(Math.random()*Ico_matrix_curr_j.size());
                            Integer dead_Ico = Ico_matrix_curr_j.get(muIco_index);
                            Ico_matrix_curr_j.remove(muIco_index);

                            icoMatrix.setDeath(dead_Ico, t_curr);

                            j++;
                        }
                        break;
                    case 20:  // death

                        minNo = Math.min(Z_curr_j, num);

                        Z_curr_j -= minNo;
                        break;
                    case 21:  //waning immunity

                        minNo = Math.min(Z_curr_j, num);

                        Z_curr_j -= minNo;
                        X_curr_j += minNo;
                        break;
                    case 22:  //migration

                        minNo = Math.min(Yco_curr_i, num);

                        Yco_curr_j += minNo;
                        Yco_curr_i -= minNo;

                        //cumultive incidence
                        cumIco_j += minNo;

                        j = 0;
                        while(j<minNo){


                            int migIco_index = (int)Math.floor(Math.random()*Ico_matrix_curr_i.size());
                            Integer mig_Ico = Ico_matrix_curr_i.get(migIco_index);


                            Ico_matrix_curr_i.remove(migIco_index);
                            //so the migrant will appear like a new infection
                            Ico_matrix_curr_j.add(Ico_matrix_length);


                            //update the entry of the migrant Ico_matrix when it was in the origin patch ;

                            //icoMatrixDeath.set(mig_Ico, t_curr);
                            icoMatrix.setDeath(mig_Ico, t_curr);

                            //this needs more thought, while migrant Ico is in another patch it will be the same coinfected individual....
                            //do we want to treat that coinfected individual differently (as a new infection), and have a different ID?
                            //It makes sense that a migrant will appear like a new infection in the new patch, but with same parents as before.
                            //And that infected individual in the original patch should be treated as 'recovered'


                            //add new entry of migrant in the new patch
                            patch = 2;

                            icoMatrix.logBirth(t_curr);//icoMatrix.getBirth(mig_Ico));
                            icoMatrix.logDeath(Double.NEGATIVE_INFINITY);
                            icoMatrix.logParent1(icoMatrix.getParent1(mig_Ico));
                            icoMatrix.logParent2(icoMatrix.getParent2(mig_Ico));
                            //icoMatrix.logFitness((icoMatrix.getFitness(mig_Ico)));
                            icoMatrix.logPatch(patch);

                            j++;
                            Ico_matrix_length++;
                        }
                    case 23:
                        // emigration of Is_i

                        minNo = Math.min(Y_curr_i, num);

                        Y_curr_j += minNo;
                        Y_curr_i -= minNo;

                        //cumulative incidence
                        cumIs_j += minNo;
                        j = 0;
                        while(j<minNo){

                            int migIs_index = (int)Math.floor(Math.random()*I_matrix_curr_i.size());
                            Integer mig_Is = I_matrix_curr_i.get(migIs_index);

                            I_matrix_curr_i.remove(migIs_index);
                            I_matrix_curr_j.add(I_matrix_length);


                            //update Ico_matrix;
                            reassortant = iMatrix.getReassortant(mig_Is);
                            if(reassortant==1) {

                                //parent_co = (Integer)iMatrixParent_co.get(mig_Is);
                                parent_co = iMatrix.getParentCo(mig_Is);
                                if(parent_co>=0) {

                                    if(Yr_primary_curr_i>0) {
                                        Yr_primary_curr_i--;
                                        Yr_primary_curr_j++;
                                        cumIrp_j++;
                                    }

                                }

                                Yr_curr_i--;
                                Yr_curr_j++;
                                cumIrs_j++;


                            }

                            death = t_curr;

                            //iMatrixDeath.set(mig_Is, death);
                            iMatrix.setDeath(mig_Is, death);

                            patch = 2;

                            iMatrix.logId(I_matrix_length);
                            iMatrix.logBirth(t_curr);//iMatrix.getBirth(mig_Is));    //i think this is wrong - it should be t_curr?
                            iMatrix.logDeath(Double.NEGATIVE_INFINITY);
                            iMatrix.logParent(mig_Is);//iMatrix.getParent(mig_Is));   //shouldn't this be mig_is, rather than its parent?
                            iMatrix.logParentCo((int) Double.NEGATIVE_INFINITY);
                            iMatrix.logReassortant(iMatrix.getReassortant(mig_Is));
                            iMatrix.logFitness(iMatrix.getFitness(mig_Is));
                            iMatrix.logSegFitness(iMatrix.getSegFitness(mig_Is));
                            iMatrix.logPatch(patch);

                            I_matrix_length++;

                            j++;
                        }
                        break;
                    case 24:

                        minNo = Math.min(X_curr_i, num);

                        X_curr_j += minNo;
                        X_curr_i -= minNo;

                        break;
                    case 25:

                        minNo = Math.min(Z_curr_i, num);

                        Z_curr_j += minNo;
                        Z_curr_i -= minNo;

                        break;
                    case 26:

                        minNo = Math.min(Yco_curr_j, num);

                        Yco_curr_i += minNo;
                        Yco_curr_j -= minNo;

                        //cumulative incidence
                        cumIco_i += minNo;


                        j = 0;
                        while(j<minNo){


                            int migIco_index = (int)Math.floor(Math.random()*Ico_matrix_curr_j.size());
                            Integer mig_Ico = Ico_matrix_curr_j.get(migIco_index);

                            Ico_matrix_curr_j.remove(migIco_index);
                            Ico_matrix_curr_i.add(Ico_matrix_length);


                            //update Ico_matrix;

                            icoMatrix.setDeath(mig_Ico, t_curr);
                            patch = 1;

                            //although the coinfected will appear like a new infection in the new patch, but it was actually generated at time currentCoinfectedHistory[1]

                            icoMatrix.logBirth(t_curr);//icoMatrix.getBirth(mig_Ico));
                            icoMatrix.logDeath(Double.NEGATIVE_INFINITY);
                            icoMatrix.logParent1(icoMatrix.getParent1(mig_Ico));
                            icoMatrix.logParent2(icoMatrix.getParent2(mig_Ico));
                            //icoMatrix.logFitness((icoMatrix.getFitness(mig_Ico)));
                            icoMatrix.logPatch(patch);

//                            if(writeOutput) {
//                                try {
//                                    writer2.write(Ico_matrix_length + "," + icoMatrix.getBirth(Ico_matrix_length) + "," + icoMatrix.getDeath(Ico_matrix_length) + "," + icoMatrix.getParent1(Ico_matrix_length) + "," +
//                                            icoMatrix.getParent2(Ico_matrix_length)+ "," + icoMatrix.getPatch(Ico_matrix_length)+ "\n");
//                                    writer2.flush();
//                                } catch (IOException e) {
//                                    e.printStackTrace();  //To change body of catch statement use File | Settings | File Templates.
//                                }
//                            }

                            j++;
                            Ico_matrix_length++;
                        }
                        break;
                    case 27:

                        //emigration of Is_j

                        minNo = Math.min(Y_curr_j, num);

                        Y_curr_i += minNo;
                        Y_curr_j -= minNo;

                        //cumulative incidence
                        cumIs_i += minNo;

                        j = 0;
                        while(j<minNo){


                            int migIs_index = (int)Math.floor(Math.random()*I_matrix_curr_j.size());
                            Integer mig_Is = I_matrix_curr_j.get(migIs_index);

                            I_matrix_curr_j.remove(migIs_index);
                            I_matrix_curr_i.add(I_matrix_length);

                            //update Ico_matrix;
                            reassortant = iMatrix.getReassortant(mig_Is);
                            if(reassortant==1) {
                                parent_co = iMatrix.getParentCo(mig_Is);
                                if(parent_co>=0) {

                                    if(Yr_primary_curr_j>0) {
                                        Yr_primary_curr_j--;
                                        Yr_primary_curr_i++;
                                        cumIrp_i++;
                                    }

                                }

                                Yr_curr_j--;
                                Yr_curr_i++;
                                cumIrs_i++;

                            }
                            death = t_curr;

                            iMatrix.setDeath(mig_Is, death);

                            patch = 1;

                            iMatrix.logId(I_matrix_length);
                            iMatrix.logBirth(t_curr);//iMatrix.getBirth(mig_Is));
                            iMatrix.logDeath(Double.NEGATIVE_INFINITY);
                            iMatrix.logParent(mig_Is);//iMatrix.getParent(mig_Is));
                            iMatrix.logParentCo((int)Double.NEGATIVE_INFINITY);
                            iMatrix.logReassortant(iMatrix.getReassortant(mig_Is));
                            iMatrix.logFitness(iMatrix.getFitness(mig_Is));
                            iMatrix.logSegFitness(iMatrix.getSegFitness(mig_Is));
                            iMatrix.logPatch(patch);

                         j++;
                            I_matrix_length++;
                        }
                        break;
                    case 28:

                        minNo = Math.min(X_curr_j, num);

                        X_curr_i += minNo;
                        X_curr_j -= minNo;

                        break;
                    case 29:

                        minNo = Math.min(Z_curr_j, num);

                        Z_curr_i += minNo;
                        Z_curr_j -= minNo;

                        break;
                    case 30:
                        //antigenic evolution in patch 1

                        j = 0;

                        while(j<num){

                            int index = (int)Math.floor(Math.random()*I_matrix_curr_i.size());
                            Integer evolvedIs = I_matrix_curr_i.get(index);
                            double currentFitness = iMatrix.getFitness(evolvedIs);

                            double newFitness = currentFitness + Math.log10((1+(params.dfe)));

                            iMatrix.setFitness(evolvedIs, newFitness);
                            iMatrix.setSegFitness(evolvedIs, newFitness);

                            j++;
                        }
                        break;
                    case 31:
                        //antigenic evolution in patch 2
                        j = 0;

                        while(j<num){

                            int index = (int)Math.floor(Math.random()*I_matrix_curr_j.size());
                            Integer evolvedIs = I_matrix_curr_j.get(index);
                            double currentFitness = iMatrix.getFitness(evolvedIs);
                            //double newFitness =currentFitness*(1+params.dfe.nextDouble());
                            double newFitness = currentFitness + Math.log10((1+(params.dfe)));

                            iMatrix.setFitness(evolvedIs, newFitness);
                            iMatrix.setSegFitness(evolvedIs, newFitness);

                            j++;
                        }
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




    private void initializeI_matrix(int Y_curr_i, int Y_curr_j) {
        for(int i=0;i<Y_curr_i;i++) {

            iMatrix.logId(i);
            iMatrix.logBirth(Double.NEGATIVE_INFINITY);
            iMatrix.logDeath(Double.NEGATIVE_INFINITY);
            iMatrix.logParent((int)Double.NEGATIVE_INFINITY);
            iMatrix.logParentCo((int)Double.NEGATIVE_INFINITY);
            iMatrix.logReassortant(0);
            iMatrix.logFitness(0.0);
            iMatrix.logSegFitness(0.0);
            iMatrix.logPatch(1);

//            if(writeOutput){
//                try {
//                    writer1.write(i + "," + iMatrix.getBirth(i) + "," + iMatrix.getDeath(i) + "," + iMatrix.getParent(i) + "," +
//                            iMatrix.getParentCo(i)+ "," + iMatrix.getReassortant(i)+ "," + iMatrix.getFitness(i)+ "," + iMatrix.getPatch(i)+ "\n");
//                } catch (IOException e) {
//                    e.printStackTrace();  //To change body of catch statement use File | Settings | File Templates.
//                }
//            }

        }

        for(int i=0;i<Y_curr_j;i++) {

            iMatrix.logId(Y_curr_i+i);
            iMatrix.logBirth(Double.NEGATIVE_INFINITY);
            iMatrix.logDeath(Double.NEGATIVE_INFINITY);
            iMatrix.logParent((int)Double.NEGATIVE_INFINITY);
            iMatrix.logParentCo((int)Double.NEGATIVE_INFINITY);
            iMatrix.logReassortant(0);
            iMatrix.logFitness(0.0);
            iMatrix.logSegFitness(0.0);
            iMatrix.logPatch(2);

//            if(writeOutput) {
//                try {
//                    writer1.write((i+Y_curr_i) + "," + iMatrix.getBirth(i+Y_curr_i) + "," + iMatrix.getDeath(i+Y_curr_i) + "," + iMatrix.getParent(i+Y_curr_i) + "," +
//                            iMatrix.getParentCo(i+Y_curr_i)+ "," + iMatrix.getReassortant(i+Y_curr_i)+ "," + iMatrix.getFitness(i+Y_curr_i)+ "," + iMatrix.getPatch(i+Y_curr_i)+ "\n");
//                } catch (IOException e) {
//                    e.printStackTrace();  //To change body of catch statement use File | Settings | File Templates.
//                }
//            }

        }

    }

    public void initializeIco_matrix(int Yco_curr_i, int Yco_curr_j, int I_curr) {
        for(int i=0;i<Yco_curr_i;i++) {

            icoMatrix.logBirth(Double.NEGATIVE_INFINITY);
            icoMatrix.logDeath(Double.NEGATIVE_INFINITY);
            icoMatrix.logParent1((int)Math.floor(Math.random()*I_curr));
            icoMatrix.logParent2((int)Math.floor(Math.random()*I_curr));
            //icoMatrix.logFitness(Distribution.nextDouble());
            icoMatrix.logPatch(1);

//            if(writeOutput) {
//                try {
//                    writer2.write(i + "," + icoMatrix.getBirth(i) + "," + icoMatrix.getDeath(i) + "," + icoMatrix.getParent1(i) + "," +
//                            icoMatrix.getParent2(i)+ "," + icoMatrix.getPatch(i)+ "\n");
//                } catch (IOException e) {
//                    e.printStackTrace();  //To change body of catch statement use File | Settings | File Templates.
//                }
//            }
        }
        for(int i=0;i<Yco_curr_j;i++) {

            icoMatrix.logBirth(Double.NEGATIVE_INFINITY);
            icoMatrix.logDeath(Double.NEGATIVE_INFINITY);
            icoMatrix.logParent1((int)Math.floor(Math.random()*I_curr));
            icoMatrix.logParent2((int)Math.floor(Math.random()*I_curr));
            //icoMatrix.logFitness(Distribution.nextDouble());
            icoMatrix.logPatch(2);

//            if(writeOutput) {
//                try {
//                    writer2.write((i+Yco_curr_i) + "," + icoMatrix.getBirth(i+Yco_curr_i) + "," + icoMatrix.getDeath(i+Yco_curr_i) + "," + icoMatrix.getParent1(i+Yco_curr_i) + "," +
//                            icoMatrix.getParent2(i+Yco_curr_i)+ "," + icoMatrix.getPatch(i+Yco_curr_i)+ "\n");
//                } catch (IOException e) {
//                    e.printStackTrace();  //To change body of catch statement use File | Settings | File Templates.
//                }
//            }
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
    public BitSet getCurrentInfected(Integer I_curr, Integer startingNumber) {

        BitSet I_matrix_curr = new BitSet();

        I_matrix_curr.set(startingNumber, I_curr);

        return I_matrix_curr;
    }
    public int[] getCoinfectedHistory(int birth, int death, int parent1, int parent2, int patch) {

        int[] coinfectHistory = new int[5];
        coinfectHistory[0] = birth;
        coinfectHistory[1] = death;
        coinfectHistory[2] = parent1;
        coinfectHistory[3] = parent2;
        coinfectHistory[4] = patch;
        return coinfectHistory;
    }
    private void processMatrix(List<int[]> I_matrix, List<int[]> Ico_matrix) {

        List<int[]> patch1_Imatrix = new ArrayList<int[]>();
        List<int[]> patch2_Imatrix = new ArrayList<int[]>();


        for(int i=0; i<I_matrix.size(); i++) {

            int[] infectionEntry = I_matrix.get(i);

            //is it a reassortant?
            int reassortant = infectionEntry[4];


            if(reassortant == 1) {

                int patch = infectionEntry[5];

                if(patch == 1) {

                    patch1_Imatrix.add(infectionEntry);

                }
                else{

                    patch2_Imatrix.add(infectionEntry);
                }

            }
        }

        System.out.println(patch1_Imatrix.size());
        System.out.println(patch2_Imatrix.size());

        patch1_Imatrix.clear();
        patch2_Imatrix.clear();
        for(int[] infectionEntry:Ico_matrix) {


            int patch = infectionEntry[4];

            if(patch == 1) {

                patch1_Imatrix.add(infectionEntry);

            }
            else{

                patch2_Imatrix.add(infectionEntry);
            }

        }

        System.out.println(patch1_Imatrix.size());
        System.out.println(patch2_Imatrix.size());

        patch1_Imatrix.clear();
        patch2_Imatrix.clear();

        for(int[] i: I_matrix) {

            int patch = i[5];

            if(patch == 1) {

                patch1_Imatrix.add(i);
            }
            else{
                patch2_Imatrix.add(i);
            }
        }

        System.out.println("patch1: "+patch1_Imatrix.size());
        System.out.println("patch2: "+patch2_Imatrix.size());

        patch1_Imatrix.clear();
        patch2_Imatrix.clear();
    }

    private List<Double> calculateWeights(List<Integer> I_matrix_curr) {

        List<Double> cumulativeWeights = new ArrayList<Double>();
        List<Double> normalisedWeights = new ArrayList<Double>();

        //Calculate the cumulativeWeights
        //Double firstWeight = (Double)iMatrixFitness.get(I_matrix_curr.get(0));
        Double firstWeight = iMatrix.getFitness(I_matrix_curr.get(0));
        cumulativeWeights.add(firstWeight);
        for(int i=1; i<I_matrix_curr.size(); i++) {

            int infected = I_matrix_curr.get(i);
            //Double weight = (Double)iMatrixFitness.get(infected) + (cumulativeWeights.get(i-1));
            Double weight = iMatrix.getFitness(infected) + (cumulativeWeights.get(i-1));
            cumulativeWeights.add(weight);

        }

        //Calculate the normalised weights
        Double totalWeight = cumulativeWeights.get(cumulativeWeights.size()-1);
        for(int i=0; i<cumulativeWeights.size(); i++) {

            Double calcWeight = (cumulativeWeights.get(i))/totalWeight;
            normalisedWeights.add(calcWeight);

        }

        return normalisedWeights;

    }
    /**
     *
     * @param cumulativeWeights normalised cumulative weights for the groups in each year
     * @return chooses a group (returns a bin number) based on cumulative probability weights
     *
     */
    public int chooseBin(List<Double> cumulativeWeights) {


        //get a number between 0 & 1 froma uniform distribution
        Double x = Math.random();

        //find the bin in which x fits e.g. if x = 0.5 and cumulativeWeights = [0.2, 0.4, 0.6, 1]
        //then final i (the bin) should be 1

        int i = 0;
        int j = cumulativeWeights.size()-1;
        while (i < j && (x > cumulativeWeights.get(i))) {
            i++;
        }
        return i;

    }

    public infectionHistory getIMatrix() {

       return iMatrix;

    }

    public coinfectionHistory getIcoMatrix() {
        return icoMatrix;
    }

    public static void main(String args[]) {

        EpiParams params = new EpiParams();
        reassortmentTwoPatch_new twoPatchModel = new reassortmentTwoPatch_new();
        twoPatchModel.runSimulation(params);

    }
}

