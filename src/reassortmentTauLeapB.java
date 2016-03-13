import cern.colt.list.DoubleArrayList;
import cern.jet.random.Poisson;
import org.apache.commons.math3.stat.descriptive.DescriptiveStatistics;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.*;

/**
 * Created with IntelliJ IDEA.
 * User: jr207
 * Date: 8/28/12
 * Time: 11:03 AM
 * To change this template use File | Settings | File Templates.
 */
public class reassortmentTauLeapB implements EpiModel {

    private double tau = 1.0/4.0;
    public infectionHistory iMatrix = new infectionHistory((int)40e6);
    public coinfectionHistory icoMatrix = new coinfectionHistory();

    private FileWriter writer1 = null;
    private FileWriter writer2 = null;
    private FileWriter writer3 = null;

    private double tmrca1;
    private double diversity1;
    private double tmrca2;
    private double diversity2;
    private double tmrca;
    private double diversity;
    private double meanPopFitness;


    private boolean writeIncid = true;




    public void runSimulation(EpiParams params, int sim_no, double [][] tmrca_seg1, double [][] tmrca_seg2, double [][] fitness, double[][] fitnessv ) {


        params.print();

        System.out.println("x");
        //icoMatrix.initialize((int)1e6);
        int patch = 1;

        File outputfile = new File("one_patch_model_reassortants_epsi_"+params.epsilon_endemic+"_D_"+params.durationOfInfection+
                          "_W_"+params.waningImmunity+"_R0_"+params.R0+"_simTime_"+params.simulationTime+"yrs_p_"+params.p+"_epistasis_"+params.e_b+"_simNo_"+(sim_no+1)+".csv");


        File cumIncid = new File("cumulative_Incidence_antigenicMu_epsi_"+params.epsilon_endemic+"_D_"+params.durationOfInfection+
                          "_W_"+params.waningImmunity+"_R0_"+params.R0+"_simTime_"+params.simulationTime+"yrs_p_"+params.p+"_epistasis_"+params.e_b+"_simNo_"+(sim_no+1)+".csv");

        //File outputfile2 = new File("tmrca_and_fitness_through_time_p_"+params.p+"_psi_"+params.psi+"_neutralMu_"+params.neutralMu+"_sn_"+params.s_n+"_epistasis_"+params.e_b+"_simNo_"+sim_no+".csv");
        try {
            writer1 = new FileWriter(outputfile);
            //writer3 = new FileWriter(outputfile2);

        } catch (IOException e) {
            e.printStackTrace();  //To change body of catch statement use File | Settings | File Templates.
        }



        if(writeIncid){
            try {
                writer2 = new FileWriter(cumIncid);
            } catch (IOException e) {
                e.printStackTrace();  //To change body of catch statement use File | Settings | File Templates.
            }

        }

        icoMatrix.initialize((int)8e6);

        double maxTime = params.simulationEndTime;

        int N = params.N;
        //initial numbers;
        int X_curr = params.S_init;
        int Y_curr = params.Is_init;
        int Yco_curr = params.Ico_init;
        int Z_curr = params.R_init;
        int Yr_curr = 0;
        int Yr_primary_curr = 0;

        DoubleArrayList rates = new DoubleArrayList();
        double t_curr = 0;
        double t_next;

        double beta_s;
        double beta_co;


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
        double evolve; // antigenic evolution
        double evolve_large;
        double weak_evolve;
        double neg_evolve;
        //double lethal_evolve;



        List<Integer> I_matrix_curr = getCurrentInfected(Y_curr, 0);  //current infections in patch i
        List<Integer> Ico_matrix_curr = getCurrentInfected(Yco_curr, 0); //current coinfections in patch i


        int I_matrix_length = Y_curr;
        int Ico_matrix_length = Yco_curr;

        initializeI_matrix(Y_curr);
        initializeIco_matrix(Yco_curr, Y_curr);

        //cumulative incidence;

        int cumIs = 0;
        int cumIco = 0;
        int cumIrp = 0;
        int cumIrs = 0;

        int p = 0;
        double sampleTime = 0;
        DescriptiveStatistics stats1 = new DescriptiveStatistics();
        DescriptiveStatistics stats2 = new DescriptiveStatistics();
        DescriptiveStatistics stats3 = new DescriptiveStatistics();

        tmrca_seg1[sim_no] = new double[61];
        tmrca_seg2[sim_no] = new double[61];
        fitness[sim_no] = new double[61];
        fitnessv[sim_no] = new double[61];
        int ind = 0;

        while (t_curr <= maxTime) {
            double popFitness = getPopFitness(I_matrix_curr);
            meanPopFitness = popFitness;
            double popFitnessVariance = getPopFitnessVar(popFitness, I_matrix_curr);


            //System.out.println(Math.pow(10,Double.NEGATIVE_INFINITY));

            if(t_curr==sampleTime) {
                updateDiversity(I_matrix_curr);


                System.out.println(t_curr + "," + (60-((t_curr-tmrca1)/365.25)) + "," + diversity1 / 365.25+"," + (60-((t_curr-tmrca2)/365.25)) + "," + diversity2 / 365.25+","+popFitness+","+popFitnessVariance+","+(tmrca1-tmrca2)/365.25);
                sampleTime +=365.25;


                tmrca_seg1[sim_no][ind] = (t_curr-tmrca1)/365.25;
                tmrca_seg2[sim_no][ind] = (t_curr-tmrca2)/365.25;
                fitness[sim_no][ind] = popFitness;
                fitnessv[sim_no][ind] = popFitnessVariance;
                ind++;


                if(t_curr/365.25 > 30) {
                    stats1.addValue(tmrca1/365.25);
                    stats2.addValue(tmrca2/365.25);
                    stats3.addValue((tmrca1-tmrca2)/365.25);

                }
            }


            //System.out.println(t_curr + "," + Math.pow(10, popFitness) + ", " + Math.pow(10, popFitnessVariance));
            try {


                writer1.write(t_curr+","+X_curr+","+Y_curr+","+Yco_curr+","+Yr_primary_curr+","+Yr_curr+","+Z_curr+","+popFitness+","+ popFitnessVariance+"\n");
                writer1.flush();
            } catch (IOException e) {
                e.printStackTrace();  //To change body of catch statement use File | Settings | File Templates.
            }



            rates.removeAll(rates);

            beta_s =  params.beta_s*(1.0+params.epsilon*Math.sin((2.0*Math.PI/365.0)*t_curr));
            beta_co = beta_s;

            muN = params.mu*N;
            BSIs = beta_s*X_curr*Y_curr/N;
            BSIco = beta_s*X_curr*Yco_curr/N;
            muS = params.mu*X_curr;
            BIsIs = params.psi*beta_co*Y_curr*Y_curr/N;
            nuIs = params.nu_s*Y_curr;
            nuIco = params.nu_co*Yco_curr;
            muIs = params.mu*Y_curr;
            muIco = params.mu*Yco_curr;
            muR = params.mu*Z_curr;
            gammaR = params.gamma*Z_curr;
            evolve = params.antigenicMu*Y_curr; //ignoring evolution in coinfected since such small numbers comparatively
            neg_evolve = params.nonAntigenicMu*Y_curr;
//            evolve_large = (params.q)*params.antigenicMu_b*Y_curr;
//            weak_evolve = params.neutralMu*Y_curr;
            //lethal_evolve = 0.001*Y_curr;

            //System.out.println("evolve_large:"+evolve_large+" q "+params.q+" p "+params.p);


            //neutral fitness effect - 1/Ne (population of infected)
            //params.s_n = 1.0/(double)Y_curr;

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
            rates.add(evolve);
            rates.add(neg_evolve);
//            rates.add(evolve_large);
//            rates.add(weak_evolve);



            t_next = tau+t_curr;
            Poisson poisson;

            int event=0;
            int noOfRates = rates.size();
            int j;

            int parent_co;
            int reassortant;

            int minNo;

            List<Double> fitnessWeights;

            if(writeIncid) {

                try {
                    writer2.write(t_curr+","+cumIs+","+cumIco+","+cumIrp+","+cumIrs+"\n");
                    writer2.flush();
                } catch (IOException e) {
                    e.printStackTrace();  //To change body of catch statement use File | Settings | File Templates.
                }

            }
            while(event<noOfRates) {
                poisson = new Poisson(tau*rates.get(event), params.randomGenerator);
                int num = poisson.nextInt();

                //System.out.println(">"+num);
                switch(event) {

                    case 0: // birth of S
                        X_curr += num;


                        break;
                    case 1: // transmission from Is

                        if(Y_curr == 0){

                            break;
                        }

                        minNo = Math.min(X_curr, num);

                        Y_curr += minNo;
                        X_curr -= minNo;

                        cumIs += minNo;
                        j = 0;

                        int I_index;
                        //Integer parent = (int)Double.NEGATIVE_INFINITY;

                        fitnessWeights = calculateWeights(I_matrix_curr);




                        while(j<minNo) {

                            //I want to choose the parent index according to the fitness weights


                            int bin = chooseParent(fitnessWeights);//fitnessWeights.indexOf(Collections.max(fitnessWeights));//(int)Math.floor(Math.random()*I_matrix_curr.size());//
                            Integer parent = I_matrix_curr.get(bin);


                            //Integer parent = topFitnesses.get(j);

                            //find out if the infected host carries a reassortant strain
                            reassortant = iMatrix.getReassortant(parent);//
                            if(reassortant == 1) {
                                Yr_curr++;
                                //cumulative Irs_i (secondary Ir)
                                cumIrs++;
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
                            iMatrix.logSeg1parent(parent);
                            iMatrix.logSeg2parent(parent);

                            I_matrix_curr.add(I_matrix_length);
                            I_matrix_length++;

                            j++;

                        }

                        break;
                    case 2: // transmission from Ico

                        if(Yco_curr==0){
                            break;
                        }

                        minNo = Math.min(X_curr, num);

                        Y_curr += minNo;
                        X_curr -= minNo;

                        //cumuative incidence
                        cumIs += minNo;


                        j = 0;

                        while(j<minNo) {

                            int Ico_index = (int)Math.floor(Math.random()*Ico_matrix_curr.size());
                            parent_co = Ico_matrix_curr.get(Ico_index);

                            int parent1 = icoMatrix.getParent1(parent_co);
                            int parent2 = icoMatrix.getParent2(parent_co);

                            if(Double.isInfinite(iMatrix.getFitness(parent1))||Double.isInfinite(iMatrix.getFitness(parent2))) {

                                while(Double.isInfinite(iMatrix.getFitness(parent1))||Double.isInfinite(iMatrix.getFitness(parent2))) {

                                    Ico_index = (int)Math.floor(Math.random()*Ico_matrix_curr.size());
                                    parent_co = Ico_matrix_curr.get(Ico_index);
                                    parent1 = icoMatrix.getParent1(parent_co);
                                    parent2 = icoMatrix.getParent2(parent_co);

                                }

                            }
                            Integer parent = (int)Double.NEGATIVE_INFINITY;

                            //how long has coinfected host been transmitting? Thus, tau = birth-t_curr and we're not using an average tau i.e 1/(Bco*S/N)

                            //double transmTime = t_curr-(Double)icoMatrixBirth.get(parentCo);  //So it is scaled in days
                            double transmTime = t_curr-icoMatrix.getBirth(parent_co);  //So it is scaled in days
//                            double nTrials = params.mutnRate*transmTime;
//                            double reassortantLost = 1.0-1.0/(double)params.viralLoad;
//                            double reassortantNotFixedInTau = Math.pow(reassortantLost,nTrials);
//                            double rho = 1.0 - reassortantNotFixedInTau;
                            double rho = 0.5;
                            double reassortantNotFixedInTau = 1-rho;

                            DoubleArrayList coinfectedTransmission = new DoubleArrayList();

                            coinfectedTransmission.add(rho);
                            coinfectedTransmission.add(reassortantNotFixedInTau/2);
                            coinfectedTransmission.add(reassortantNotFixedInTau/2);

                            double random = Math.random();
                            int virus = chooseEvent(coinfectedTransmission, random, 1.0);

                            double newFitness = 0.0;double seg1Fitness = 0.0;int seg1parent = (int)Double.NEGATIVE_INFINITY; int seg2parent = (int)Double.NEGATIVE_INFINITY;
                            if(virus>0) {

                                reassortant = 0;
                                if(virus==1){

                                    parent = icoMatrix.getParent1(parent_co);//Integer)icoMatrixParent1.get(parentCo);
                                    newFitness = iMatrix.getFitness(parent);//+errorDistribution.nextDouble();
                                    seg1Fitness = iMatrix.getSegFitness(parent);
                                    seg1parent = parent;
                                    seg2parent = parent;

                                }
                                else{
//
                                    parent = icoMatrix.getParent2(parent_co);//Integer)icoMatrixParent2.get(parentCo);
                                    newFitness = iMatrix.getFitness(parent);//+errorDistribution.nextDouble();
                                    seg1Fitness = iMatrix.getSegFitness(parent);
                                    seg1parent = parent;
                                    seg2parent = parent;
                                }

                            }
                            else{

                                reassortant = 1;

                                double parent1Fitness = iMatrix.getFitness(icoMatrix.getParent1(parent_co));
                                double parent2Fitness = iMatrix.getFitness(icoMatrix.getParent2(parent_co));

                                double seg1Fitness_p1 = iMatrix.getSegFitness(icoMatrix.getParent1(parent_co));
                                double seg1Fitness_p2 = iMatrix.getSegFitness(icoMatrix.getParent2(parent_co));

                                double seg2Fitness_p1 = parent1Fitness-seg1Fitness_p1;
                                double seg2Fitness_p2 = parent2Fitness-seg1Fitness_p2;



                                //only need to do this when segment 2 is neutrally evolving, as the viral fitness is solely governed by
                                // segment 1, so need to choose which of the coinfected parent donates segment 1 (randomly)

                                int r = (int)Math.round(Math.random());

                                //newFitness = (parent1Fitness + parent2Fitness)*0.5 + errorDistribution.nextDouble();
                                //newFitness = ((parent1Fitness + errorDistribution.nextDouble()) + (parent2Fitness + errorDistribution.nextDouble()))*0.5;

                                if(r < 1) {  // seg 1 inherited from parent 1 and seg 2 inherited from parent 2


                                    newFitness = seg1Fitness_p1 + seg2Fitness_p2;
                                    seg1Fitness = seg1Fitness_p1;
                                    seg1parent = icoMatrix.getParent1(parent_co);
                                    seg2parent = icoMatrix.getParent2(parent_co);
                                }
                                else{       // seg 1 inherited from parent 2 and seg 2 inherited from parent 1

                                    newFitness = seg1Fitness_p2 + seg2Fitness_p1;
                                    seg1Fitness = seg1Fitness_p2;
                                    seg1parent = icoMatrix.getParent2(parent_co);
                                    seg2parent = icoMatrix.getParent1(parent_co);

                                }


                                //newFitness = (parent1Fitness + parent2Fitness)*0.5;

                                Yr_primary_curr++;
                                Yr_curr++;

                                //cumulative Irp and Irs in patch i
                                cumIrp++;
                                cumIrs++;

                            }

                            if(Double.isNaN(newFitness)) {
                                System.out.println(5);
                            }
                            iMatrix.logId(I_matrix_length);
                            iMatrix.logBirth(t_curr);
                            iMatrix.logDeath(Double.NEGATIVE_INFINITY);
                            iMatrix.logParent(parent);
                            iMatrix.logParentCo(parent_co);
                            iMatrix.logReassortant(reassortant);
                            iMatrix.logFitness(newFitness);
                            iMatrix.logSegFitness(seg1Fitness);
                            iMatrix.logSeg1parent(seg1parent);
                            iMatrix.logSeg2parent(seg2parent);

                            //iMatrix.logPatch(patch);


                            I_matrix_curr.add(I_matrix_length);
                            I_matrix_length++;


                            j++;
                        }



                        break;
                    case 3: // death of S

                        minNo = Math.min(X_curr, num);
                        X_curr -= minNo;

                        break;

                    case 4: // transmission from Is to another Is - coinfection



                        minNo = Math.min(Y_curr, num);

                        Yco_curr += minNo;
                        Y_curr -= minNo;

                        //cumulative Ico incid
                        cumIco += minNo;


                        j = 0;

                        while(j<minNo){

                            int recipient_index = (int)Math.floor(Math.random()*I_matrix_curr.size());
                            Integer recipientParent = I_matrix_curr.get(recipient_index);
                            I_matrix_curr.remove(recipientParent);

                            //donor virus/host

                            int donor_index = (int)Math.floor(Math.random()*I_matrix_curr.size());
                            int donorParent = I_matrix_curr.get(donor_index);

                            reassortant = iMatrix.getReassortant(recipientParent);
                            //should I have update the I_matrix for the recipient virus - will appear as if it is dead
                            if(reassortant == 1) {

                                parent_co = iMatrix.getParentCo(recipientParent);
                                if(parent_co >= 0) {

                                    if(Yr_primary_curr>0){
                                        Yr_primary_curr--;
                                    }

                                }

                                Yr_curr--;
                            }

                            icoMatrix.logBirth(t_curr);
                            icoMatrix.logDeath(Double.NEGATIVE_INFINITY);
                            icoMatrix.logParent1(donorParent);
                            icoMatrix.logParent2(recipientParent);
                            //icoMatrix.logFitness(fitnessDistribution.nextDouble());
                            //icoMatrix.logPatch(patch);

                            Ico_matrix_curr.add(Ico_matrix_length);
                            Ico_matrix_length++;

                            j++;

                        }




                        break;
                    case 5: // recovery of Is

                        minNo = Math.min(Y_curr, num);
                        Y_curr -= minNo;
                        Z_curr += minNo;

                        j = 0;
                        while(j<minNo){

                            int nuIs_index = (int)Math.floor(Math.random()*I_matrix_curr.size());
                            Integer recovered_Is = I_matrix_curr.get(nuIs_index);
                            I_matrix_curr.remove(nuIs_index);

                            reassortant = iMatrix.getReassortant(recovered_Is);
                            if(reassortant==1) {

                                parent_co = iMatrix.getParentCo(recovered_Is);
                                if(parent_co >= 0) {
                                    //check if primary reassortant
                                    if(Yr_primary_curr>0) {
                                        Yr_primary_curr--;
                                    }

                                }

                                Yr_curr--;
                            }

                            iMatrix.setDeath(recovered_Is, t_curr);

                            j++;
                        }

                        break;
                    case 6: // recovery of Ico

                        minNo = Math.min(Yco_curr, num);

                        Yco_curr -= minNo;
                        Z_curr += minNo;

                        j = 0;
                        while(j<minNo) {

                            int nuIco_index = (int)Math.floor(Math.random()*Ico_matrix_curr.size());
                            Integer recovered_Ico = Ico_matrix_curr.get(nuIco_index);
                            Ico_matrix_curr.remove(nuIco_index);

                            icoMatrix.setDeath(recovered_Ico, t_curr);

                            j++;
                        }

                        break;
                    case 7: // death of Is

                        minNo = Math.min(Y_curr, num);
                        Y_curr -= minNo;
                        j = 0;
                        while(j<minNo) {


                            int muIs_index = (int)Math.floor(Math.random()*I_matrix_curr.size());
                            Integer dead_Is = I_matrix_curr.get(muIs_index);
                            I_matrix_curr.remove(muIs_index);

                            reassortant = iMatrix.getReassortant(dead_Is);
                            if(reassortant==1) {

                                parent_co = iMatrix.getParentCo(dead_Is);
                                if(parent_co>=0) {

                                    if(Yr_primary_curr>0) {
                                        Yr_primary_curr--;
                                    }

                                }

                                Yr_curr--;
                            }

                            iMatrix.setDeath(dead_Is, t_curr);


                            j++;
                        }

                        break;
                    case 8: // death of Ico

                        minNo = Math.min(Yco_curr, num);

                        Yco_curr -= minNo;

                        j = 0;
                        while(j<minNo) {


                            int muIco_index = (int)Math.floor(Math.random()*Ico_matrix_curr.size());
                            Integer dead_Ico = Ico_matrix_curr.get(muIco_index);
                            Ico_matrix_curr.remove(muIco_index);

                            icoMatrix.setDeath(dead_Ico, t_curr);

                            j++;
                        }
                        break;
                    case 9: // death of R

                        minNo = Math.min(Z_curr, num);
                        Z_curr -= minNo;

                        break;
                    case 10: // waning immunity

                        minNo = Math.min(Z_curr, num);

                        Z_curr -= minNo;
                        X_curr += minNo;

                        break;

                    case 11:

                        j = 0;

                        // seg 1 evolution
                        while(j<num){

                            int index = (int)Math.floor(Math.random()*I_matrix_curr.size());
                            Integer evolvedIs = I_matrix_curr.get(index);
                            double currentFitness = iMatrix.getFitness(evolvedIs);
                            //double newFitness = currentFitness *((1+(params.dfe.nextDouble())));

                            if(Double.isInfinite(currentFitness)) {
                                while(Double.isInfinite(currentFitness)) {

                                    index = (int)Math.floor(Math.random()*I_matrix_curr.size());
                                    evolvedIs = I_matrix_curr.get(index);
                                    currentFitness = iMatrix.getFitness(evolvedIs);
                                }
                            }
                            //double newFitness = currentFitness + Math.pow(Math.log10(1 + params.s_b), Math.exp(params.e_b));


//                            double newFitness = currentFitness + (Math.log10(1 + params.benDFE.nextDouble()));
//
//                            if(Double.isNaN(newFitness)) {
//                                System.out.println(1);
//                            }

                            //double newFitness = currentFitness + Math.log10((1+(params.expDFE.nextDouble())));

                            double e_c = Math.random();

                            double newFitness = currentFitness;
                            if(e_c <= params.p_ben1) {

                                newFitness += Math.log10(1+params.benDFE.nextDouble());
                            }
                            else{
                                newFitness +=  Math.log10(1-params.delDFE.nextDouble()); //+ Math.log10(1+epistasis_2.nextDouble());
                            }



                            iMatrix.setFitness(evolvedIs, newFitness);
                            iMatrix.setSegFitness(evolvedIs, newFitness);   // keeps track of seg1 fitness, when neutral these total and seg fitnesses are the same

                            j++;
                        }
                        break;

                    case 12:

                        // non-antigenic mutation on seg2

                        j = 0;
                        Collections.shuffle(I_matrix_curr);

                        while(j<num){

                            //int index = (int)Math.floor(Math.random()*I_matrix_curr.size());
                            int index = 0;
                            Integer evolvedIs = I_matrix_curr.get(j);
                            double currentFitness = iMatrix.getFitness(evolvedIs);

                            if(Double.isInfinite(currentFitness)) {
                                while(Double.isInfinite(currentFitness)) {

                                    index = (int)Math.floor(Math.random()*I_matrix_curr.size());
                                    evolvedIs = I_matrix_curr.get(index);
                                    currentFitness = iMatrix.getFitness(evolvedIs);
                                }
                            }

                            double e_c = Math.random();

                            double newFitness = currentFitness;
                            if(e_c <= params.p_ben2) {

                                newFitness += Math.log10(1+params.benDFE.nextDouble());
                            }
                            else{
                                newFitness +=  Math.log10(1-params.delDFE.nextDouble()); //+ Math.log10(1+epistasis_2.nextDouble());
                            }



                            iMatrix.setFitness(evolvedIs, newFitness);


                            j++;
                        }
                        break;


                }
                event++;
            }
            t_curr = t_next;
        }

        iMatrix.patch = new ArrayList<Integer>(Collections.nCopies(iMatrix.birth.size(),patch));

        double meanTMRCA1 = stats1.getMean();
        double meanTMRCA2 = stats2.getMean();
        double meanDiff = stats3.getMean();

        double std = stats1.getStandardDeviation();
        double min = stats1.getMin();
        double max = stats1.getMax();
        double lowerq = stats1.getPercentile(25);
        double upperq = stats1.getPercentile(75);
        double median = stats1.getPercentile(50);
        double std2 = stats2.getStandardDeviation();
        double min2 = stats2.getMin();
        double max2 = stats2.getMax();
        double lowerq2 = stats2.getPercentile(25);
        double upperq2 = stats2.getPercentile(75);
        double median2 = stats2.getPercentile(50);

        double std3 = stats3.getStandardDeviation();
        double min3 = stats3.getMin();
        double max3 = stats3.getMax();
        double lowerq3 = stats3.getPercentile(25);
        double upperq3 = stats3.getPercentile(75);
        double median3 = stats3.getPercentile(50);

        System.out.println("mean1: "+ meanTMRCA1);
        System.out.println("std: "+ std);
        System.out.println("min: "+ min+", max: "+max);
        System.out.println("LQ: "+ lowerq+", UQ: "+upperq);
        System.out.println();
        System.out.println("mean2: "+ meanTMRCA2);
        System.out.println("std: "+ std2);
        System.out.println("min: "+ min2+", max: "+max2);
        System.out.println("LQ: "+ lowerq2+", UQ: "+upperq2);
        System.out.println();
        System.out.println("mean3: "+ meanDiff);
        System.out.println("std: "+ std3);
        System.out.println("min: "+ min3+", max: "+max3);
        System.out.println("LQ: "+ lowerq3+", UQ: "+upperq3);


        //System.out.println("><"+I_matrix_length);
    }

    public double getPopFitness(List<Integer> I_curr) {

        Double totalFitness = 0.0;
        int count = 0;
        for(Integer i: I_curr) {


            if(iMatrix.getFitness(i) > 0) {
                totalFitness += iMatrix.getFitness((i));
                count++;
            }
        }


        return (totalFitness/count);


    }

    public double getPopFitnessVar(double popFitness, List<Integer> I_curr) {

        Double totalSquaredDeviations = 0.0;

        int count = 0;
        for(Integer i: I_curr) {

            double fitness = iMatrix.getFitness(i);
            if(fitness != Double.NEGATIVE_INFINITY && fitness > 0) {
                totalSquaredDeviations += Math.pow((fitness - popFitness), 2);
            }
            else{
                count++;
            }
        }

        return (totalSquaredDeviations)/(I_curr.size()-count);

    }

    private int commonAncestor(int i1, int i2) {

        int commonAnc = (int)Double.NEGATIVE_INFINITY;
        int lineageA = i1;
        int lineageB = i2;
        Set<Integer> ancestry = new HashSet<Integer>();
        while(true) {

            if(iMatrix.getParent(lineageA)!=(int)Double.NEGATIVE_INFINITY) {

                lineageA = iMatrix.getParent(lineageA);
                if (!ancestry.add(lineageA)) {
                    commonAnc = lineageA;
                    break;
                }
            }
            if(iMatrix.getParent(lineageB)!=(int)Double.NEGATIVE_INFINITY) {

                lineageB = iMatrix.getParent(lineageB);
                if (!ancestry.add(lineageB)) {
                    commonAnc = lineageB;
                    break;
                }
            }
            if (iMatrix.getParent(lineageA) == (int)Double.NEGATIVE_INFINITY && iMatrix.getParent(lineageB) == (int)Double.NEGATIVE_INFINITY) {
                break;
            }

        }
        return commonAnc;
    }

    private int commonAncestor(int i1, int i2, int seg) {

        int commonAnc = (int)Double.NEGATIVE_INFINITY;
        int lineageA = i1;
        int lineageB = i2;
        Set<Integer> ancestry = new HashSet<Integer>();
        if(seg==0) {
            while (true) {

                if (iMatrix.getSeg1parent(lineageA) != (int) Double.NEGATIVE_INFINITY) {

                    lineageA = iMatrix.getSeg1parent(lineageA);
                    if (!ancestry.add(lineageA)) {
                        commonAnc = lineageA;
                        break;
                    }
                }
                if (iMatrix.getSeg1parent(lineageB) != (int) Double.NEGATIVE_INFINITY) {

                    lineageB = iMatrix.getSeg1parent(lineageB);
                    if (!ancestry.add(lineageB)) {
                        commonAnc = lineageB;
                        break;
                    }
                }
                if (iMatrix.getSeg1parent(lineageA) == (int) Double.NEGATIVE_INFINITY && iMatrix.getSeg1parent(lineageB) == (int) Double.NEGATIVE_INFINITY) {
                    break;
                }

            }
        }
        else{
            while (true) {

                if (iMatrix.getSeg2parent(lineageA) != (int) Double.NEGATIVE_INFINITY) {

                    lineageA = iMatrix.getSeg2parent(lineageA);
                    if (!ancestry.add(lineageA)) {
                        commonAnc = lineageA;
                        break;
                    }
                }
                if (iMatrix.getSeg2parent(lineageB) != (int) Double.NEGATIVE_INFINITY) {

                    lineageB = iMatrix.getSeg2parent(lineageB);
                    if (!ancestry.add(lineageB)) {
                        commonAnc = lineageB;
                        break;
                    }
                }
                if (iMatrix.getSeg2parent(lineageA) == (int) Double.NEGATIVE_INFINITY && iMatrix.getSeg2parent(lineageB) == (int) Double.NEGATIVE_INFINITY) {
                    break;
                }

            }
        }
        return commonAnc;
    }

    private double distance(int i1, int i2) {


        int ancestor = commonAncestor(i1,i2);
        if(ancestor!= (int)Double.NEGATIVE_INFINITY) {

            double distA = iMatrix.getBirth(i1) - iMatrix.getBirth(ancestor);
            double distB = iMatrix.getBirth(i2) - iMatrix.getBirth(ancestor);
            return distA + distB;
        }
        else{
            return 0;

        }
    }

    private double distance(int i1, int i2, int seg) {


        int ancestor = commonAncestor(i1,i2,seg);
        if(ancestor!= (int)Double.NEGATIVE_INFINITY) {

            double distA = iMatrix.getBirth(i1) - iMatrix.getBirth(ancestor);
            double distB = iMatrix.getBirth(i2) - iMatrix.getBirth(ancestor);
            return distA + distB;
        }
        else{
            return 0;

        }
    }


    public void updateDiversity(List<Integer> I_matrix_curr) {
        diversity1 = 0.0;
        tmrca1 = 0.0;
        diversity2 = 0.0;
        tmrca2 = 0.0;
        int sampleCount1 = 0;
        int sampleCount2 = 0;



        List<Integer> indices = new ArrayList<Integer>();

        for(int i= 0; i < I_matrix_curr.size(); i++) {

            indices.add(i);
        }

        Collections.shuffle(indices, new Random((long)(Math.floor(Math.random()*10000))));
        int samplingDepth = I_matrix_curr.size()/100;

        for (int i = 0; i < 50; i++) {


            int j = samplingDepth*i;
            //System.out.println(j);
            Integer vA = I_matrix_curr.get(i);
            //int r = (int)Math.floor(Math.random()*I_matrix_curr.size());
            Integer vB = I_matrix_curr.get(indices.get(i+1));
            if (vA != (int)Double.NEGATIVE_INFINITY && vB != (int)Double.NEGATIVE_INFINITY) {
                double dist = distance(vA, vB, 0);
                diversity1 += dist;
                if (dist > tmrca1) {
                    tmrca1 = dist;
                }
                sampleCount1 += 1;
            }

            if (vA != (int)Double.NEGATIVE_INFINITY && vB != (int)Double.NEGATIVE_INFINITY) {
                double dist = distance(vA, vB, 1);
                diversity2 += dist;
                if (dist > tmrca2) {
                    tmrca2 = dist;
                }
                sampleCount2 += 1;
            }
        }
        if (sampleCount1 > 0) {
            diversity1 /= (double) sampleCount1;
        }
        if (sampleCount2 > 0) {
            diversity2 /= (double) sampleCount2;

        }
        tmrca1 /= 2.0;
        tmrca2 /= 2.0;

    }


    public void updateDiversity(List<Integer> I_matrix_curr, int seg) {
        diversity = 0.0;
        tmrca = 0.0;
        int sampleCount = 0;


        int samplingDepth = I_matrix_curr.size()/200;

        for (int i = 0; i < 200; i++) {


            int j = samplingDepth*i;
            //System.out.println(j);
            Integer vA = I_matrix_curr.get(j);
            Integer vB = I_matrix_curr.get(j+1);
            if (vA != (int)Double.NEGATIVE_INFINITY && vB != (int)Double.NEGATIVE_INFINITY) {
                double dist = distance(vA, vB, seg);
                diversity += dist;
                if (dist > tmrca) {
                    tmrca = dist;
                }
                sampleCount += 1;
            }
        }
        if (sampleCount > 0) {
            diversity /= (double) sampleCount;
        }
        tmrca /= 2.0;
    }


    @Override
    public void runSimulation(EpiParams params) {

    }

    @Override
    public double getRatesSum(DoubleArrayList rates) {
        return 0;  //To change body of implemented methods use File | Settings | File Templates.
    }

    public List<Integer> getCurrentInfected(int Y_curr, int startingNumber) {

        List<Integer> I_matrix_curr = new ArrayList<Integer>();
        int i = startingNumber;
        while(i<Y_curr){

            I_matrix_curr.add(i);
            i++;
        }

        return I_matrix_curr;

    }

    private void initializeI_matrix(int Y_curr) {


        for(int i=0;i<Y_curr;i++) {

            iMatrix.logId(i);
            iMatrix.logBirth(Double.NEGATIVE_INFINITY);
            iMatrix.logDeath(Double.NEGATIVE_INFINITY);
            iMatrix.logParent((int)Double.NEGATIVE_INFINITY);
            iMatrix.logParentCo((int)Double.NEGATIVE_INFINITY);
            iMatrix.logReassortant(0);
            iMatrix.logSeg1parent((int)Double.NEGATIVE_INFINITY);
            iMatrix.logSeg2parent((int)Double.NEGATIVE_INFINITY);


            double fitness = Math.log10(1.2);
            iMatrix.logFitness(fitness);
            iMatrix.logSegFitness(fitness);



        }

//        for(int i=0;i<100;i++) {
//
//            iMatrix.logId(i);
//            iMatrix.logBirth(Double.NEGATIVE_INFINITY);
//            iMatrix.logDeath(Double.NEGATIVE_INFINITY);
//            iMatrix.logParent((int)Double.NEGATIVE_INFINITY);
//            iMatrix.logParentCo((int)Double.NEGATIVE_INFINITY);
//            iMatrix.logReassortant(0);
//
//            double fitness = Math.log10(1.2);
//            iMatrix.logFitness(fitness);
//            iMatrix.logSegFitness(fitness);
// }
    }

    public void initializeIco_matrix(int Yco_curr, int I_curr) {
        for(int i=0;i<Yco_curr;i++) {

            icoMatrix.logBirth(Double.NEGATIVE_INFINITY);
            icoMatrix.logDeath(Double.NEGATIVE_INFINITY);
            icoMatrix.logParent1((int)Math.floor(Math.random()*I_curr));
            icoMatrix.logParent2((int)Math.floor(Math.random()*I_curr));


        }
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
    private void processMatrix(List<int[]> I_matrix, List<int[]> Ico_matrix) {

        List<int[]> patch_Imatrix = new ArrayList<int[]>();


        for(int i=0; i<I_matrix.size(); i++) {

            int[] infectionEntry = I_matrix.get(i);

            //is it a reassortant?
            int reassortant = infectionEntry[4];
            //System.out.println(reassortant);


            if(reassortant == 1 && infectionEntry[3]>=0) {

                patch_Imatrix.add(infectionEntry);

            }
        }

        System.out.println(">"+patch_Imatrix.size());


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


    private List<Double> calculateWeights(List<Integer> I_matrix_curr) {

        List<Double> cumulativeWeights = new ArrayList<Double>();
        List<Double> fitnessValues = new ArrayList<Double>();
        List<Double> relativeFitnessValues = new ArrayList<Double>();

        Collections.sort(I_matrix_curr, comparator);


        for(Integer i: I_matrix_curr) {

            fitnessValues.add(iMatrix.getFitness(i));
        }

        Double minFit = Collections.min(fitnessValues);
        Double maxFit = Collections.max(fitnessValues);


        for(Integer i: I_matrix_curr) {

            Double fitness = iMatrix.getFitness(i);


            Double relfitness = fitness/maxFit;//(Math.pow(10, fitness-maxFit));//Math.log10(Math.pow(10, fitness/maxFit));//(fitness/maxFit)-1.0);


            relativeFitnessValues.add(relfitness);
        }
        int count = Collections.frequency(relativeFitnessValues,0.0);

        return relativeFitnessValues;
        //return relativeFitnessValues;

    }


    private List<Double> calculateWeights(List<Integer> I_matrix_curr, double meanPopFitness) {

        final List<Double> relativeFitnessValues = new ArrayList<Double>();
        List<Double> cumulativeWeights = new ArrayList<Double>();
        List<Double> normalisedWeights = new ArrayList<Double>();



        //Calculate the cumulativeWeights
        // Double firstWeight = (Double)iMatrix.Fitness.get(I_matrix_curr.get(0));
//
        Double firstWeight = iMatrix.getFitness(I_matrix_curr.get(0));
//        cumulativeWeights.add(firstWeight/meanPopFitness);
        for (int i = 0; i < I_matrix_curr.size(); i++) {

            int infected = I_matrix_curr.get(i);
            //Double weight = (Double)iMatrixFitness.get(infected) + (cumulativeWeights.get(i-1));
            Double relativeFitness = (iMatrix.getFitness(infected));
            relativeFitnessValues.add(relativeFitness);

        }

        Double maxRelFit = Collections.max(relativeFitnessValues);
        if(relativeFitnessValues.get(0) < 0.0) {
            cumulativeWeights.add((Math.abs(relativeFitnessValues.get(0)) / (Math.pow(maxRelFit,2)))) ;

        }
        else {
            cumulativeWeights.add(relativeFitnessValues.get(0) / maxRelFit);
        }
        for (int i = 1; i < relativeFitnessValues.size(); i++) {


            if(relativeFitnessValues.get(i) < 0.0) {
                cumulativeWeights.add((Math.abs(relativeFitnessValues.get(i)) / (Math.pow(maxRelFit, 2))) + cumulativeWeights.get(i - 1));
                //System.out.println(maxRelFit+","+relativeFitnessValues.get(i)+","+(Math.abs(relativeFitnessValues.get(i)) / (Math.pow(maxRelFit,2))) + cumulativeWeights.get(i - 1));
            }
            else {
                cumulativeWeights.add((relativeFitnessValues.get(i) / maxRelFit) + cumulativeWeights.get(i - 1));
            }
        }


        //       }

        //Calculate the normalised weights
        double totalWeight = cumulativeWeights.get(cumulativeWeights.size() - 1);
        for (Double cumulativeWeight : cumulativeWeights) {

            double calcWeight = cumulativeWeight / (float) totalWeight;
            normalisedWeights.add(calcWeight);
        }


        //Collections.sort(I_matrix_curr, comparator);

        double max = Collections.max(normalisedWeights);
        return normalisedWeights;
    }


    private  Comparator<Integer>  comparator = new Comparator<Integer>() {
        @Override
        public int compare(Integer o1, Integer o2) {


            double diff = iMatrix.getFitness(o1) - iMatrix.getFitness(o2);
            if(diff == 0) {

                return 0;
            }
            else if(diff > 0) {
                return 1;
            }
            else{
                return -1;
            }
        }
    };




    /**
     *
     * @param cumulativeWeights normalised cumulative weights for the groups in each year
     * @return chooses a group (returns a bin number) based on cumulative probability weights
     *
     */
    public int chooseBin(List<Double> cumulativeWeights) {


        //get a number between 0 & 1 from a uniform distribution
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

    public int chooseParent(List<Double> relativeFitnesses) {



        Integer parent_index = (int)Math.floor(Math.random() * relativeFitnesses.size());

        Double parent_fitness = relativeFitnesses.get(parent_index);

        if(parent_fitness == 1.0) {

            return parent_index;
        }
        else{

            Double x = Math.random();
            while(parent_fitness < x ) {
                //System.out.println(parent_fitness+","+x);

                parent_index = (int)Math.floor(Math.random()*relativeFitnesses.size());
                parent_fitness = relativeFitnesses.get(parent_index);
                if(parent_fitness == 1.0) {
                    break;
                }
                x = Math.random();
            }

            Double real_fitness = iMatrix.getFitness(parent_index);
            return parent_index;

        }

    }
    public infectionHistory getIMatrix() {
        return iMatrix;
    }

    public coinfectionHistory getIcoMatrix() {
        return icoMatrix;
    }

    public static void main(String [] args) {

        EpiParams params = new EpiParams();




//        params.antigenicMu = Double.parseDouble(args[0]);
//        params.s_b = Double.parseDouble(args[1]);

        params.psi = Double.parseDouble(args[0]);

        boolean getTrees = Boolean.parseBoolean(args[1]);
        int no_of_sims = Integer.parseInt(args[2]);

        params.antigenicMu = Double.parseDouble(args[3]);
        params.s_ben = Double.parseDouble(args[4]);
        params.nonAntigenicMu = Double.parseDouble(args[5]);
        params.s_del = Double.parseDouble(args[6]);
        params.p_ben1 = Double.parseDouble(args[7]);
        params.p_ben2 = Double.parseDouble(args[8]);


        double[][] tmrca_seg1 = new double[no_of_sims][61];
        double[][] tmrca_seg2 = new double[no_of_sims][61];
        double[][] fitness = new double[no_of_sims][61];
        double[][] fitnessv = new double[no_of_sims][61];


        for(int i=0; i < no_of_sims; i++) {

            reassortmentTauLeapB model = new reassortmentTauLeapB();
            model.runSimulation(params, i, tmrca_seg1, tmrca_seg2, fitness, fitnessv);

            SimulateTree tree = new SimulateTree();
//            tree.n_lineages = 100;
//
            if(getTrees) {
                tree.sampleStartTime = 30 * 365.25;
                tree.sampleEndTime = params.simulationEndTime;
                tree.n_lineages = 300;
                //System.out.println(tree.samplingSchemeForOnePatchModel);
                tree.getTransmissionTrees(2, model, (i+1));
            }


        }

        File tmrca_seg1_output = new File("tmrca1_through_time_antigenicMu_"+params.antigenicMu+"_nonAntigenicMu_"+params.nonAntigenicMu+"_sb_"
                                        +params.s_ben+"_sd_"+params.s_del+"_propB1_"+params.p_ben1+"_propB2_"+params.p_ben2+".csv");
        File tmrca_seg2_output = new File("tmrca2_and_through_time_antigenicMu_"+params.antigenicMu+"_nonAntigenicMu_"+params.nonAntigenicMu+"_sb_"
                                        +params.s_ben+"_sd_"+params.s_del+"_propB1_"+params.p_ben1+"_propB2_"+params.p_ben2+".csv");
        File fitness_output = new File("fitness_through_time_antigenicMu_"+params.antigenicMu+"_nonAntigenicMu_"+params.nonAntigenicMu+"_sb_"
                                        +params.s_ben+"_sd_"+params.s_del+"_propB1_"+params.p_ben1+"_propB2_"+params.p_ben2+".csv");

        try {
            BufferedWriter writer1 = new BufferedWriter(new FileWriter(tmrca_seg1_output));
            BufferedWriter writer2 = new BufferedWriter(new FileWriter(tmrca_seg2_output));
            BufferedWriter writer3 = new BufferedWriter(new FileWriter(fitness_output));
            //BufferedWriter writer4 = new BufferedWriter(new FileWriter(fitnessv_output));

            for(int i=0; i < 61; i++) {

                String tmrca1_time = "";
                String tmrca1_raw = "";
                String tmrca2_time = "";
                String tmrca2_raw = "";
                String fitness_string = "";
                String fitnessv_string = "";


                for (int j = 0; j < no_of_sims; j++) {

                    tmrca1_raw += tmrca_seg1[j][i] + ",";
                    tmrca2_raw += tmrca_seg2[j][i] + ",";
                    tmrca1_time += (60 - tmrca_seg1[j][i]) + ",";
                    tmrca2_time += (60 - tmrca_seg2[j][i]) + ",";
                    fitness_string += fitness[j][i] + ",";
                    fitnessv_string += fitnessv[j][i] + ",";


                }

                writer1.write(i + "," + tmrca1_time+tmrca1_raw + "\n");
                writer2.write(i + "," + tmrca2_time+tmrca2_raw + "\n");
                writer3.write(i + "," + fitness_string+fitnessv_string + "\n");

                writer1.flush();
                writer2.flush();
                writer3.flush();
            }
        } catch (IOException e) {
            e.printStackTrace();
        }



    }

}

