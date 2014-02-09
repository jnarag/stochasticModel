import cern.colt.list.DoubleArrayList;
import cern.jet.random.Normal;
import cern.jet.random.Poisson;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collections;
import java.util.List;

/**
 * Created with IntelliJ IDEA.
 * User: jr207
 * Date: 8/28/12
 * Time: 11:03 AM
 * To change this template use File | Settings | File Templates.
 */
public class reassortmentTauLeap implements EpiModel {

    private double tau = 1.0/4.0;
    public infectionHistory iMatrix = new infectionHistory((int)40e6);
    public coinfectionHistory icoMatrix = new coinfectionHistory();

    private FileWriter writer1 = null;
    private FileWriter writer2 = null;


    private boolean writeIncid = true;

    public void runSimulation(EpiParams params) {

        params.print();

        System.out.println("x");
        //icoMatrix.initialize((int)1e6);
        int patch = 1;

        File outputfile = new File("one_patch_model_reassortants_antigenicMu_"+params.antigenMu+"_s_"+params.dfe+"_psi_"+params.psi+"_rho_"+params.mutnRate+"_"+"_epsi_"+params.epsilon_endemic+
                "_D_"+params.durationOfInfection+"_W_"+params.waningImmunity+"_R0_"+params.R0+"_simTime_"+params.simulationTime+"yrs.csv");


        File cumIncid = new File("cumulative_Incidence_antigenicMu_"+params.antigenMu+"_s_"+params.dfe+"psi_"+params.psi+"_rho_"+params.mutnRate+"_"+"_epsi_"+params.epsilon_endemic+
                "_D_"+params.durationOfInfection+"_W_"+params.waningImmunity+"_R0_"+params.R0+"_simTime_"+params.simulationTime+"yrs.csv");

        try {
            writer1 = new FileWriter(outputfile);
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


        List<Integer> I_matrix_curr = getCurrentInfected(Y_curr, 0);  //current infections in patch i
        List<Integer> Ico_matrix_curr = getCurrentInfected(Yco_curr, 0); //current coinfections in patch i

        //Fitness parameters;
        // Normal fitnessDistribution = new Normal(6.0, 2.0, params.randomGenerator);
        Normal errorDistribution = new Normal(0.0, 2.0, params.randomGenerator);

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
        while (t_curr <= maxTime) {

            System.out.println(t_curr+","+X_curr+","+Y_curr+","+Yco_curr+","+Yr_primary_curr+","+Yr_curr+","+Z_curr);

            //printPopFitness(I_matrix_curr);

            try {

                double popFitness = getPopFitness(I_matrix_curr);
                double popFitnessVariance = getPopFitnessVar(popFitness, I_matrix_curr);

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
            evolve = params.antigenMu*Y_curr; //ignoring evolution in coinfected since such small numbers comparatively

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
            //rates.add(evolve);
            //rates.add(0.0); // seg 2 evolution is same as seg 1

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

                            int bin = chooseBin(fitnessWeights);
                            Integer parent = I_matrix_curr.get(bin);

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

                            Integer parent = (int)Double.NEGATIVE_INFINITY;

                            //how long has coinfected host been transmitting? Thus, tau = birth-t_curr and we're not using an average tau i.e 1/(Bco*S/N)

                            //double transmTime = t_curr-(Double)icoMatrixBirth.get(parentCo);  //So it is scaled in days
                            double transmTime = t_curr-icoMatrix.getBirth(parent_co);  //So it is scaled in days
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

                                    parent = icoMatrix.getParent1(parent_co);//Integer)icoMatrixParent1.get(parentCo);
                                    newFitness = iMatrix.getFitness(parent);//+errorDistribution.nextDouble();
                                    seg1Fitness = iMatrix.getSegFitness(parent);
                                }
                                else{
//
                                    parent = icoMatrix.getParent2(parent_co);//Integer)icoMatrixParent2.get(parentCo);
                                    newFitness = iMatrix.getFitness(parent);//+errorDistribution.nextDouble();
                                    seg1Fitness = iMatrix.getSegFitness(parent);
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
                                }
                                else{       // seg 2 inherited from parent 2 and seg 2 inherited from parent 1

                                    newFitness = seg1Fitness_p2 + seg2Fitness_p1;
                                    seg1Fitness = seg1Fitness_p2;

                                }


                                //newFitness = (parent1Fitness + parent2Fitness)*0.5;

                                Yr_primary_curr++;
                                Yr_curr++;

                                //cumulative Irp and Irs in patch i
                                cumIrp++;
                                cumIrs++;

                            }

                            iMatrix.logId(I_matrix_length);
                            iMatrix.logBirth(t_curr);
                            iMatrix.logDeath(Double.NEGATIVE_INFINITY);
                            iMatrix.logParent(parent);
                            iMatrix.logParentCo(parent_co);
                            iMatrix.logReassortant(reassortant);
                            iMatrix.logFitness(newFitness);
                            iMatrix.logSegFitness(seg1Fitness);
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
                            double newFitness = currentFitness + Math.log10((1+(params.dfe)));
                            //double newFitness = currentFitness + Math.log10((1+(params.expDFE.nextDouble())));


                            iMatrix.setFitness(evolvedIs, newFitness);
                            iMatrix.setSegFitness(evolvedIs, newFitness);   // keeps track of seg1 fitness, when neutral these total and seg fitnesses are the same
//                            double currentFitness = iMatrix.getFitness(evolvedIs);
//                            double newFitness = currentFitness + Math.log10((1+(params.dfe.nextDouble())));
//                            iMatrix.setFitness(evolvedIs, newFitness);

                            j++;
                        }
                        break;

                    case 12:

                        // seg 2 evolution currently the same mutation rate and dfe as seg 1

                        j = 0;

                        while(j<num){

                            int index = (int)Math.floor(Math.random()*I_matrix_curr.size());
                            Integer evolvedIs = I_matrix_curr.get(index);
                            double currentFitness = iMatrix.getFitness(evolvedIs);
                            //double newFitness = currentFitness *((1+(params.dfe.nextDouble())));
                            double newFitness = currentFitness + Math.log10((1+(params.dfe)));


                            //iMatrix.setSegFitness(evolvedIs, iMatrix.getFitness(evolvedIs));
                            //segFitness doesn't improve...
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

        //System.out.println("><"+I_matrix_length);
    }

    public double getPopFitness(List<Integer> I_curr) {

        Double totalFitness = 0.0;
        for(Integer i: I_curr) {

            totalFitness += iMatrix.getFitness(i);
        }

        //System.out.println("Pop fitness: "+ totalFitness/I_curr.size());

        return totalFitness/I_curr.size();


    }

    public double getPopFitnessVar(double popFitness, List<Integer> I_curr) {

        Double totalSquaredDeviations = 0.0;

        for(Integer i: I_curr) {

            totalSquaredDeviations += Math.pow((iMatrix.getFitness(i)-popFitness), 2);

        }

        return totalSquaredDeviations/I_curr.size();

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
            iMatrix.logFitness(0.0);
            iMatrix.logSegFitness(0.0);


        }
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
//
//        System.out.println(Ico_matrix.size());
//
//
//
//        System.out.println(I_matrix.size());

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

    public static void main(String [] args) {

        EpiParams params = new EpiParams();
        reassortmentTauLeap model = new reassortmentTauLeap();
        model.runSimulation(params);

        infectionHistory iMatrix = model.getIMatrix();

//        List<int[]> I_matrix = model.I_matrix;
//        List<int[]> Ico_matrix = model.Ico_matrix;

        //model.processMatrix(I_matrix,Ico_matrix);

        System.out.println("Imat size: "+iMatrix.birth.size()+ " : "+iMatrix.death.size()+ " : "+iMatrix.parent.size()+" : "+iMatrix.parentCo.size()+" : "+iMatrix.reassortant.size());


    }

}


