import cern.colt.list.DoubleArrayList;
import cern.jet.random.Exponential;
import cern.jet.random.Poisson;
import org.apache.commons.math3.stat.descriptive.DescriptiveStatistics;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.math.RoundingMode;
import java.text.DecimalFormat;
import java.util.*;

/**
 * Created by jayna on 23/07/2016.
 */
public class moranReassortment implements EpiModel {


    private double tau = 1.0/4.0;
    public infectionHistory iMatrix = new infectionHistory((int)40e6);
    public coinfectionHistory icoMatrix = new coinfectionHistory();

    private FileWriter writer1 = null;
    private FileWriter writer2 = null;
    private FileWriter writer3 = null;
    private FileWriter writer4 = null;

    private double tmrca1;
    private double diversity1;
    private double tmrca2;
    private double diversity2;
    private double tmrca;
    private double diversity;
    private double meanPopFitness;

    private boolean writeIncid = true;







    public void runSimulation(EpiParams params, int sim_no, double [][] tmrca_seg1, double [][] tmrca_seg2, double [][] fitness, double[][] fitnessv, double [][] fitness_1, double[][] fitnessv_1, double [][] fitness_2, double[][] fitnessv_2) {




        Map<Double, Map<String, Double>> fitDistTime = new HashMap<Double, Map<String, Double>>();

//      writer1//        if(writeIncid){
//            try {
//                writer2 = new FileWriter(cumIncid);
//            } catch (IOException e) {
//                e.printStackTrace();  //To change body of catch statement use File | Settings | File Templates.
//            }
//
//        }


        icoMatrix.initialize((int) 8e6);

        //Exponential benDFE_2 = new Exponential(1/0.3, params.randomGenerator);
        //Exponential delDFE_2 = new Exponential(1/sd1, params.randomGenerator);

        //System.out.println("sd1 "+sd1);

        double maxTime = params.simulationEndTime;
        DoubleArrayList rates = new DoubleArrayList();

        int Y_curr = params.Is_init;
        int Yco_curr = params.Ico_init;
        int Yr_curr = 0;
        int Yr_primary_curr = 0;

        double t_curr = 0;
        double t_next;

        double BSIs;
        double BIsIs;
        double BSIco;


        List<Integer> I_matrix_curr = getCurrentInfected(Y_curr, 0);  //current infections in patch i
        List<Integer> Ico_matrix_curr = getCurrentInfected(Yco_curr, 0); //current coinfections in patch i
        int I_matrix_length = Y_curr;
        int Ico_matrix_length = Yco_curr;

        initializeI_matrix(Y_curr);
        initializeIco_matrix(Yco_curr, Y_curr);

        double sampleTime = 0;
        double sampleTime1 = 0;

        DescriptiveStatistics stats1 = new DescriptiveStatistics();
        DescriptiveStatistics stats2 = new DescriptiveStatistics();
        DescriptiveStatistics stats3 = new DescriptiveStatistics();

        DescriptiveStatistics var = new DescriptiveStatistics();

        tmrca_seg1[sim_no] = new double[61];
        tmrca_seg2[sim_no] = new double[61];

        fitness[sim_no] = new double[61];
        fitnessv[sim_no] = new double[61];

        fitness_1[sim_no] = new double[61];
        fitnessv_1[sim_no] = new double[61];

        fitness_2[sim_no] = new double[61];
        fitnessv_2[sim_no] = new double[61];
        int ind = 0;


        //cumulative incidence;

        int cumIs = 0;
        int cumIco = 0;
        int cumIrp = 0;
        int cumIrs = 0;


        while (t_curr <= maxTime) {

            double popFitness = getPopFitness(I_matrix_curr);
            meanPopFitness = popFitness;
            double popFitnessVariance = getPopFitnessVar(popFitness, I_matrix_curr);

            if (t_curr == sampleTime) {
                updateDiversity(I_matrix_curr);

                //System.out.println(t_curr + "," + Y_curr + "," + Yco_curr + "," + Yr_primary_curr + "," + (t_curr - tmrca1) / 365 + "," + diversity1 / 365 + "," + (tmrca1 - tmrca2) / 365 + "," + popFitness+", "+popFitnessVariance);

                //System.out.println(t_curr + "," + (60-((t_curr-tmrca1)/365)) + "," + diversity1 / 365+"," + (60-((t_curr-tmrca2)/365)) + "," + diversity2 / 365+","+popFitness+","+popFitnessVariance+","+(tmrca1-tmrca2)/365+","+Y_curr+","+Yco_curr);
                sampleTime += 365;



//                if (t_curr / 365 > 40) {
//
//                    System.out.println(t_curr / 365);
                    double seg1_Fitness_mean = getSegPopFitness(I_matrix_curr, 0);
                    double seg2_Fitness_mean = getSegPopFitness(I_matrix_curr, 1);
                    double seg1_Fitness_var = getPopFitnessVar(seg1_Fitness_mean, I_matrix_curr, 0);
                    double seg2_Fitness_var = getPopFitnessVar(seg2_Fitness_mean, I_matrix_curr, 1);

                    tmrca_seg1[sim_no][ind] = (t_curr - tmrca1) / 365;
                    tmrca_seg2[sim_no][ind] = (t_curr - tmrca2) / 365;

                    fitness[sim_no][ind] = popFitness;
                    fitnessv[sim_no][ind] = popFitnessVariance;

                    fitness_1[sim_no][ind] = seg1_Fitness_mean;
                    fitnessv_1[sim_no][ind] = seg1_Fitness_var;

                    fitness_2[sim_no][ind] = seg2_Fitness_mean;
                    fitnessv_2[sim_no][ind] = seg2_Fitness_var;
                    ind++;

                    stats1.addValue(tmrca1 / 365);
                    stats2.addValue(tmrca2 / 365);
                    stats3.addValue((tmrca1 - tmrca2) / 365);
                    var.addValue(popFitnessVariance);

                    //System.out.println("mean "+seg1_Fitness_mean+", "+seg2_Fitness_mean+", "+popFitness);
                    //System.out.println("var "+seg1_Fitness_var+", "+seg2_Fitness_var+", "+popFitnessVariance);


//                }
            }

            if (t_curr == sampleTime1) {

                sampleTime1 += 182.5;

                List<Double> fitnessValues = new ArrayList<Double>();

                for (Integer i : I_matrix_curr) {

                    fitnessValues.add(iMatrix.getFitness(i));
                }

                Collections.sort(fitnessValues);

                DecimalFormat df = new DecimalFormat("#.####");
                df.setRoundingMode(RoundingMode.CEILING);


                Double min = Collections.min(fitnessValues);
                Double max = Collections.max(fitnessValues);

//                String firstBin = df.format(min-0.1);
//                String lastBin = df.format(max+0.1);

                Double firstBin = (Math.round((min - 0.1000000000) * 10000.0) / 10000.0);
                Double lastBin = (Math.round((max + 0.1000000000) * 10000.0) / 10000.0);

                List<Double> fitBins = new ArrayList<Double>();

                Double bin = firstBin + 0.1;//Double.parseDouble(firstBin)+0.1;

                fitBins.add(firstBin);
                while (bin < lastBin) {

                    fitBins.add(bin);
                    bin += 0.1;
                }

                Map<String, Double> fitBinsFreq = new HashMap<String, Double>();

                double total_freq = 0;
                double total_count = 0;
                for (Double b : fitBins) {

                    Double fitbin = b;

                    Double lower_bin = fitbin;
                    Double upper_bin = (Math.round((fitbin + 0.1000000000) * 10000.0) / 10000.0);
                    //System.out.println(lower_bin+"-"+upper_bin);

                    double count = 0;
                    for (int f = 0; f < fitnessValues.size(); f++) {

                        Double ff = fitnessValues.get(f);

                        if (ff <= upper_bin && ff > lower_bin) {

                            count++;

                        }
                    }
                    total_count += count;

                    total_freq += count / I_matrix_curr.size();
                    fitBinsFreq.put(df.format(b), count / I_matrix_curr.size());

                }

                //if(t_curr >= 70*365) {
                //    System.out.println(total_freq + "," + total_count + "," + min + "," + max + "," + firstBin + "," + lastBin);
                //}
                fitDistTime.put(t_curr, fitBinsFreq);


            }
//
//            try {
//
//                writer1.write(t_curr+","+Y_curr+","+Yco_curr+","+Yr_primary_curr+","+Yr_curr+","+popFitness+","+ popFitnessVariance+","+diversity1/365+","+tmrca1/365+","+tmrca2/365+"\n");
//                writer1.flush();
//            } catch (IOException e) {
//                e.printStackTrace();  //To change body of catch statement use File | Settings | File Templates.
//            }

            rates.removeAll(rates);

            BSIs = 0.25 * (Y_curr + Yco_curr);//params.beta_s*params.S_init* (Y_curr+Yco_curr)/params.N; //birth/death rate

            BIsIs = params.psi * (0.25 / params.relativeRate) * (Y_curr);//params.psi*params.beta_co* Y_curr*Y_curr/params.N; //coinfection rate

            rates.add(BSIs);
            rates.add(BIsIs);

            t_next = tau + t_curr;
            Poisson poisson;

            int event = 0;
            int noOfRates = rates.size();
            int j;

            Integer parent_co;
            int reassortant;

            int minNo;

            List<Double> fitnessWeights;

            //make this simpler
//            if(writeIncid) {
//
//                try {
//                    writer2.write(t_curr+","+cumIs+","+cumIco+","+cumIrp+","+cumIrs+"\n");
//                    writer2.flush();
//                } catch (IOException e) {
//                    e.printStackTrace();  //To change body of catch statement use File | Settings | File Templates.
//                }
//
//            }

            while (event < noOfRates) {

                poisson = new Poisson(tau * rates.get(event), params.randomGenerator);
                int num = poisson.nextInt();

                //System.out.println(t_curr+","+ num+","+BSIs+","+BIsIs+","+BSIco);


                switch (event) {

                    case 0: // birth Is

                        if (Y_curr == 0) {

                            break;
                        }

                        minNo = Math.min(params.S_init, num);

//                        Y_curr += minNo;
//                        X_curr -= minNo;

                        cumIs += minNo;
                        j = 0;

                        fitnessWeights = calculateWeights(I_matrix_curr);


                        double fraction_coinfected = ((double) Yco_curr) / ((double) (Yco_curr + Y_curr));

                        while (j < minNo) {

                            //I want to choose the parent index according to the fitness weights
                            Y_curr++;

                            double Is_or_ico = Math.random();

                            if(Yco_curr == 0.0) {

                                Is_or_ico = 1.0;
                            }

                            if (Is_or_ico > fraction_coinfected) {

                                int bin = chooseParent(fitnessWeights);//fitnessWeights.indexOf(Collections.max(fitnessWeights));//(int)Math.floor(Math.random()*I_matrix_curr.size());//
                                Integer parent = I_matrix_curr.get(bin);


                                //Integer parent = topFitnesses.get(j);

                                //find out if the infected host carries a reassortant strain
                                reassortant = iMatrix.getReassortant(parent);//
                                if (reassortant == 1) {
                                    Yr_curr++;
                                    //cumulative Irs_i (secondary Ir)
                                    cumIrs++;
                                }

                                Double newFitness = iMatrix.getFitness(parent);

                                //fitness change occurs at transmission

                                Double parent_birth_time = iMatrix.getBirth(parent);
                                if (Double.isInfinite(parent_birth_time)) {
                                    parent_birth_time = 0.0;
                                }
                                //Double circulating_time = t_curr-parent_birth_time;

                                Integer mut = Poisson.staticNextInt(params.rawMutationRate);


                                double prob_b1 = params.p_b_seg1 * 0.5;
                                double prob_b2 = params.p_b_seg2 * 0.5;
                                double prob_d1 = params.p_d_seg1 * 0.5;
                                double prob_d2 = params.p_d_seg2 * 0.5;


                                Integer n_b_seg1 = 0;
                                Integer n_d_seg1 = 0;
                                Integer n_b_seg2 = 0;
                                Integer n_d_seg2 = 0;

                                for (int i = 0; i < mut; i++) {

                                    double r = Math.random();

                                    if (r <= prob_b1) {
                                        n_b_seg1++;
                                    } else if (r > prob_b1 && r <= (prob_b1 + prob_d1)) {
                                        n_d_seg1++;
                                    } else if (r <= (prob_b1 + prob_d1 + prob_b2) && r > (prob_b1 + prob_d1)) {
                                        n_b_seg2++;
                                    } else if(r > (prob_b1 + prob_d1 + prob_b2) && r <= (prob_b1 + prob_d1 + prob_b2 + prob_d2)) {
                                        n_d_seg2++;
                                    }

                                }


                                //System.out.println((mut1+mut2)+","+(b_mut_seg1+b_mut_seg2+d_mut_seg1+d_mut_seg2));

                                double seg1fitness_new = 0.0;
                                double seg2fitness_new = 0.0;

                                for (int s1 = 0; s1 < n_b_seg1; s1++) {

                                    double sb = params.benDFE.nextDouble();

                                    if (sb < 0 || sb > 1) {

                                        while (sb < 0 || sb > 1) {

                                            sb = params.benDFE.nextDouble();
                                            //System.out.println("sb "+sb);


                                        }
                                    }

                                    seg1fitness_new += Math.log10(1 + sb);

                                }
                                for (int s1 = 0; s1 < n_d_seg1; s1++) {

                                    double sd = params.delDFE.nextDouble();

                                    if (sd < 0 || sd > 1) {

                                        while (sd < 0 || sd > 1) {

                                            sd = params.delDFE.nextDouble();

                                        }
                                    }

                                    seg1fitness_new += Math.log10(1 - sd);
                                }

                                for (int s2 = 0; s2 < n_b_seg2; s2++) {

                                    double sb = params.benDFE.nextDouble();

                                    if (sb < 0 || sb > 1) {

                                        while (sb < 0 || sb > 1) {

                                            sb = params.benDFE.nextDouble();
                                            //System.out.println("sb "+sb);


                                        }
                                    }

                                    seg2fitness_new += Math.log10(1 + sb);
                                }

                                for (int s2 = 0; s2 < n_d_seg2; s2++) {

                                    double sd = params.delDFE.nextDouble();//params.delDFE.nextDouble();

                                    if (sd < 0 || sd > 1) {

                                        while (sd < 0 || sd > 1) {

                                            sd = params.delDFE.nextDouble();//params.delDFE.nextDouble();

                                        }
                                    }

                                    seg2fitness_new += Math.log10(1 - sd);

                                }


                                newFitness += (seg1fitness_new + seg2fitness_new);


                                //log infection in the corresponding vectors;
                                iMatrix.logId(I_matrix_length);
                                iMatrix.logBirth(t_curr);
                                iMatrix.logDeath(Double.NEGATIVE_INFINITY);
                                iMatrix.logParent(parent);
                                iMatrix.logParentCo((int) Double.NEGATIVE_INFINITY);
                                iMatrix.logReassortant(reassortant);
                                iMatrix.logFitness(newFitness);
                                iMatrix.logSegFitness(iMatrix.getSegFitness(parent) + seg1fitness_new);
                                iMatrix.logSeg1parent(parent);
                                iMatrix.logSeg2parent(parent);

                                I_matrix_curr.add(I_matrix_length);
                                I_matrix_length++;

                            } else {

                                //coinfected transmission

                                try {

                                    int Ico_index = (int) Math.floor(Math.random() * Ico_matrix_curr.size());
                                    parent_co = Ico_matrix_curr.get(Ico_index);
                                }
                                catch(Exception e){

                                    System.out.println(Is_or_ico+","+fraction_coinfected);
                                    System.out.println(Ico_matrix_curr.size()+","+Yco_curr+","+Y_curr);

                                    throw new RuntimeException(e);
                                }

                                int Ico_index = (int) Math.floor(Math.random() * Ico_matrix_curr.size());
                                parent_co = Ico_matrix_curr.get(Ico_index);

                                int parent1 = icoMatrix.getParent1(parent_co);
                                int parent2 = icoMatrix.getParent2(parent_co);

                                if (Double.isInfinite(iMatrix.getFitness(parent1)) || Double.isInfinite(iMatrix.getFitness(parent2))) {

                                    while (Double.isInfinite(iMatrix.getFitness(parent1)) || Double.isInfinite(iMatrix.getFitness(parent2))) {

                                        Ico_index = (int) Math.floor(Math.random() * Ico_matrix_curr.size());
                                        parent_co = Ico_matrix_curr.get(Ico_index);
                                        parent1 = icoMatrix.getParent1(parent_co);
                                        parent2 = icoMatrix.getParent2(parent_co);

                                    }

                                }
                                Integer parent = (int) Double.NEGATIVE_INFINITY;

                                //how long has coinfected host been transmitting? Thus, tau = birth-t_curr and we're not using an average tau i.e 1/(Bco*S/N)
                                double rho = 0.5;
                                double reassortantNotFixedInTau = 1 - rho;

                                DoubleArrayList coinfectedTransmission = new DoubleArrayList();

                                coinfectedTransmission.add(rho);
                                coinfectedTransmission.add(reassortantNotFixedInTau / 2);
                                coinfectedTransmission.add(reassortantNotFixedInTau / 2);

                                double random = Math.random();
                                int virus = chooseEvent(coinfectedTransmission, random, 1.0);

                                double newFitness = 0.0;
                                double seg1Fitness = 0.0;
                                int seg1parent = (int) Double.NEGATIVE_INFINITY;
                                int seg2parent = (int) Double.NEGATIVE_INFINITY;
                                if (virus > 0) {

                                    reassortant = 0;
                                    if (virus == 1) {

                                        parent = icoMatrix.getParent1(parent_co);//Integer)icoMatrixParent1.get(parentCo);
                                        newFitness = iMatrix.getFitness(parent);//+errorDistribution.nextDouble();
                                        seg1Fitness = iMatrix.getSegFitness(parent);
                                        seg1parent = parent;
                                        seg2parent = parent;

                                    } else {
//
                                        parent = icoMatrix.getParent2(parent_co);//Integer)icoMatrixParent2.get(parentCo);
                                        newFitness = iMatrix.getFitness(parent);//+errorDistribution.nextDouble();
                                        seg1Fitness = iMatrix.getSegFitness(parent);
                                        seg1parent = parent;
                                        seg2parent = parent;
                                    }

                                } else {

                                    reassortant = 1;

                                    double parent1Fitness = iMatrix.getFitness(icoMatrix.getParent1(parent_co));
                                    double parent2Fitness = iMatrix.getFitness(icoMatrix.getParent2(parent_co));

                                    double seg1Fitness_p1 = iMatrix.getSegFitness(icoMatrix.getParent1(parent_co));
                                    double seg1Fitness_p2 = iMatrix.getSegFitness(icoMatrix.getParent2(parent_co));

                                    double seg2Fitness_p1 = parent1Fitness - seg1Fitness_p1;
                                    double seg2Fitness_p2 = parent2Fitness - seg1Fitness_p2;


                                    //only need to do this when segment 2 is neutrally evolving, as the viral fitness is solely governed by
                                    // segment 1, so need to choose which of the coinfected parent donates segment 1 (randomly)

                                    int r = (int) Math.round(Math.random());

                                    //newFitness = (parent1Fitness + parent2Fitness)*0.5 + errorDistribution.nextDouble();
                                    //newFitness = ((parent1Fitness + errorDistribution.nextDouble()) + (parent2Fitness + errorDistribution.nextDouble()))*0.5;

                                    if (r < 1) {  // seg 1 inherited from parent 1 and seg 2 inherited from parent 2


                                        newFitness = seg1Fitness_p1 + seg2Fitness_p2;
                                        seg1Fitness = seg1Fitness_p1;
                                        seg1parent = icoMatrix.getParent1(parent_co);
                                        seg2parent = icoMatrix.getParent2(parent_co);
                                    } else {       // seg 1 inherited from parent 2 and seg 2 inherited from parent 1

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

                                if (Double.isNaN(newFitness)) {
                                    System.out.println(5);
                                }

                                Integer mut = Poisson.staticNextInt(params.rawMutationRate);


                                double prob_b1 = params.p_b_seg1 * 0.5;
                                double prob_b2 = params.p_b_seg2 * 0.5;
                                double prob_d1 = params.p_d_seg1 * 0.5;
                                double prob_d2 = params.p_d_seg2 * 0.5;


                                Integer n_b_seg1 = 0;
                                Integer n_d_seg1 = 0;
                                Integer n_b_seg2 = 0;
                                Integer n_d_seg2 = 0;

                                for (int i = 0; i < mut; i++) {

                                    double r = Math.random();

                                    if (r <= prob_b1) {
                                        n_b_seg1++;
                                    } else if (r > prob_b1 && r <= (prob_b1 + prob_d1)) {
                                        n_d_seg1++;
                                    } else if (r <= (prob_b1 + prob_d1 + prob_b2) && r > (prob_b1 + prob_d1)) {
                                        n_b_seg2++;
                                    } else if(r > (prob_b1 + prob_d1 + prob_b2) && r <= (prob_b1 + prob_d1 + prob_b2 + prob_d2)) {
                                        n_d_seg2++;
                                    }

                                }


                                //System.out.println((mut1+mut2)+","+(b_mut_seg1+b_mut_seg2+d_mut_seg1+d_mut_seg2));

                                double seg1fitness_new = 0.0;
                                double seg2fitness_new = 0.0;

                                for (int s1 = 0; s1 < n_b_seg1; s1++) {

                                    double sb = params.benDFE.nextDouble();

                                    if (sb < 0 || sb > 1) {

                                        while (sb < 0 || sb > 1) {

                                            sb = params.benDFE.nextDouble();
                                            //System.out.println("sb "+sb);


                                        }
                                    }

                                    seg1fitness_new += (Math.log10(1 + sb));

                                }
                                for (int s1 = 0; s1 < n_d_seg1; s1++) {

                                    double sd = params.delDFE.nextDouble();

                                    if (sd < 0 || sd > 1) {

                                        while (sd < 0 || sd > 1) {

                                            sd = params.delDFE.nextDouble();

                                        }
                                    }


                                    seg1fitness_new += (Math.log10(1 - sd));
                                }

                                for (int s2 = 0; s2 < n_b_seg2; s2++) {

                                    double sb = params.benDFE.nextDouble();

                                    if (sb < 0 || sb > 1) {

                                        while (sb < 0 || sb > 1) {

                                            sb = params.benDFE.nextDouble();
                                            //System.out.println("sb "+sb);


                                        }
                                    }


                                    seg2fitness_new += (Math.log10(1 + sb));
                                }

                                for (int s2 = 0; s2 < n_d_seg2; s2++) {

                                    double sd = params.delDFE.nextDouble();//params.delDFE.nextDouble();

                                    if (sd < 0 || sd > 1) {

                                        while (sd < 0 || sd > 1) {

                                            sd = params.delDFE.nextDouble();//params.delDFE.nextDouble();

                                        }
                                    }


                                    seg2fitness_new += (Math.log10(1 - sd));

                                }


                                newFitness += (seg1fitness_new + seg2fitness_new);


                                iMatrix.logId(I_matrix_length);
                                iMatrix.logBirth(t_curr);
                                iMatrix.logDeath(Double.NEGATIVE_INFINITY);
                                iMatrix.logParent(parent);
                                iMatrix.logParentCo(parent_co);
                                iMatrix.logReassortant(reassortant);
                                iMatrix.logFitness(newFitness);
                                iMatrix.logSegFitness(seg1Fitness + seg1fitness_new);
                                iMatrix.logSeg1parent(seg1parent);
                                iMatrix.logSeg2parent(seg2parent);


                                I_matrix_curr.add(I_matrix_length);
                                I_matrix_length++;


                            }

                            //chooseDeath(t_curr, I_matrix_curr, Ico_matrix_curr, Yr_curr, Yr_primary_curr, true);

                            fraction_coinfected = (double) Yco_curr / (double) (Y_curr + Yco_curr);


                            Is_or_ico = Math.random();

                            if (Is_or_ico < fraction_coinfected) {
                                int i2 = (int) Math.floor(Math.random() * Ico_matrix_curr.size());
                                Integer dead_Ico = Ico_matrix_curr.get(i2);
                                Ico_matrix_curr.remove(dead_Ico);
                                //fitnessWeights.remove(i2);

                                Yco_curr--;
                                icoMatrix.setDeath(dead_Ico, t_curr);
                            } else {
                                int i2 = (int) Math.floor(Math.random() * I_matrix_curr.size());
                                Integer dead_Is = I_matrix_curr.get(i2);
                                I_matrix_curr.remove(dead_Is);
                                //fitnessWeights.remove(i2);

                                reassortant = iMatrix.getReassortant(dead_Is);
                                parent_co = iMatrix.getParentCo(dead_Is);

                                if (reassortant == 1) {

                                    Yr_curr--;
                                    //cumulative Irs_i (secondary Ir)

                                    if (parent_co >= 0 && Yr_primary_curr > 0) {
                                        //System.out.println(">"+Yr_primary_curr);
                                        Yr_primary_curr--;
                                        //System.out.println(Yr_primary_curr);

                                    }
                                }

                                Y_curr--;
                                iMatrix.setDeath(dead_Is, t_curr);

                            }

                            j++;

                        }

                        break;

                    case 1: // coinfection


                        minNo = Math.min(Y_curr, num);

//                        Yco_curr += minNo;
//                        Y_curr -= minNo;

                        //cumulative Ico incid
                        cumIco += minNo;


                        j = 0;

                        //fitnessWeights = calculateWeights(I_matrix_curr);


                        while (j < minNo) {

                            Y_curr--;
                            Yco_curr++;

                            int recipient_index = (int) Math.floor(Math.random() * I_matrix_curr.size());
                            Integer recipientParent = I_matrix_curr.get(recipient_index);
                            I_matrix_curr.remove(recipientParent);

                            //itnessWeights.remove(recipient_index);
                            //donor virus/host
                            fitnessWeights = calculateWeights(I_matrix_curr);

                            int donor_index = chooseParent(fitnessWeights);
                            //(int)Math.floor(Math.random()*I_matrix_curr.size()); // if choosing parent randomly
                            int donorParent = I_matrix_curr.get(donor_index);


                            reassortant = iMatrix.getReassortant(recipientParent);
                            //should I have update the I_matrix for the recipient virus - will appear as if it is dead
                            if (reassortant == 1) {

                                parent_co = iMatrix.getParentCo(recipientParent);
                                if (parent_co >= 0) {

                                    if (Yr_primary_curr > 0) {
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

                    case 2: //transmission from Ico

                        if (Yco_curr == 0) {
                            break;
                        }

                        minNo = Math.min(Yco_curr, num);

                        j = 0;
//
//                        Y_curr += minNo;


                        //cumuative incidence
                        cumIs += minNo;

                        while (j < minNo) {

                            Y_curr++;

                            int Ico_index = (int) Math.floor(Math.random() * Ico_matrix_curr.size());
                            parent_co = Ico_matrix_curr.get(Ico_index);

                            int parent1 = icoMatrix.getParent1(parent_co);
                            int parent2 = icoMatrix.getParent2(parent_co);

                            if (Double.isInfinite(iMatrix.getFitness(parent1)) || Double.isInfinite(iMatrix.getFitness(parent2))) {

                                while (Double.isInfinite(iMatrix.getFitness(parent1)) || Double.isInfinite(iMatrix.getFitness(parent2))) {

                                    Ico_index = (int) Math.floor(Math.random() * Ico_matrix_curr.size());
                                    parent_co = Ico_matrix_curr.get(Ico_index);
                                    parent1 = icoMatrix.getParent1(parent_co);
                                    parent2 = icoMatrix.getParent2(parent_co);

                                }

                            }
                            Integer parent = (int) Double.NEGATIVE_INFINITY;

                            //how long has coinfected host been transmitting? Thus, tau = birth-t_curr and we're not using an average tau i.e 1/(Bco*S/N)
                            double rho = 0.5;
                            double reassortantNotFixedInTau = 1 - rho;

                            DoubleArrayList coinfectedTransmission = new DoubleArrayList();

                            coinfectedTransmission.add(rho);
                            coinfectedTransmission.add(reassortantNotFixedInTau / 2);
                            coinfectedTransmission.add(reassortantNotFixedInTau / 2);

                            double random = Math.random();
                            int virus = chooseEvent(coinfectedTransmission, random, 1.0);

                            double newFitness = 0.0;
                            double seg1Fitness = 0.0;
                            int seg1parent = (int) Double.NEGATIVE_INFINITY;
                            int seg2parent = (int) Double.NEGATIVE_INFINITY;
                            if (virus > 0) {

                                reassortant = 0;
                                if (virus == 1) {

                                    parent = icoMatrix.getParent1(parent_co);//Integer)icoMatrixParent1.get(parentCo);
                                    newFitness = iMatrix.getFitness(parent);//+errorDistribution.nextDouble();
                                    seg1Fitness = iMatrix.getSegFitness(parent);
                                    seg1parent = parent;
                                    seg2parent = parent;

                                } else {
//
                                    parent = icoMatrix.getParent2(parent_co);//Integer)icoMatrixParent2.get(parentCo);
                                    newFitness = iMatrix.getFitness(parent);//+errorDistribution.nextDouble();
                                    seg1Fitness = iMatrix.getSegFitness(parent);
                                    seg1parent = parent;
                                    seg2parent = parent;
                                }

                            } else {

                                reassortant = 1;

                                double parent1Fitness = iMatrix.getFitness(icoMatrix.getParent1(parent_co));
                                double parent2Fitness = iMatrix.getFitness(icoMatrix.getParent2(parent_co));

                                double seg1Fitness_p1 = iMatrix.getSegFitness(icoMatrix.getParent1(parent_co));
                                double seg1Fitness_p2 = iMatrix.getSegFitness(icoMatrix.getParent2(parent_co));

                                double seg2Fitness_p1 = parent1Fitness - seg1Fitness_p1;
                                double seg2Fitness_p2 = parent2Fitness - seg1Fitness_p2;


                                //only need to do this when segment 2 is neutrally evolving, as the viral fitness is solely governed by
                                // segment 1, so need to choose which of the coinfected parent donates segment 1 (randomly)

                                int r = (int) Math.round(Math.random());

                                //newFitness = (parent1Fitness + parent2Fitness)*0.5 + errorDistribution.nextDouble();
                                //newFitness = ((parent1Fitness + errorDistribution.nextDouble()) + (parent2Fitness + errorDistribution.nextDouble()))*0.5;

                                if (r < 1) {  // seg 1 inherited from parent 1 and seg 2 inherited from parent 2


                                    newFitness = seg1Fitness_p1 + seg2Fitness_p2;
                                    seg1Fitness = seg1Fitness_p1;
                                    seg1parent = icoMatrix.getParent1(parent_co);
                                    seg2parent = icoMatrix.getParent2(parent_co);
                                } else {       // seg 1 inherited from parent 2 and seg 2 inherited from parent 1

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

                            if (Double.isNaN(newFitness)) {
                                System.out.println(5);
                            }

                            Integer mut = Poisson.staticNextInt(params.rawMutationRate);


                            double prob_b1 = params.p_b_seg1 * 0.5;
                            double prob_b2 = params.p_b_seg2 * 0.5;
                            double prob_d1 = params.p_d_seg1 * 0.5;
                            double prob_d2 = params.p_d_seg2 * 0.5;


                            Integer n_b_seg1 = 0;
                            Integer n_d_seg1 = 0;
                            Integer n_b_seg2 = 0;
                            Integer n_d_seg2 = 0;

                            for (int i = 0; i < mut; i++) {

                                double r = Math.random();

                                if (r <= prob_b1) {
                                    n_b_seg1++;
                                } else if (r > prob_b1 && r <= (prob_b1 + prob_d1)) {
                                    n_d_seg1++;
                                } else if (r <= (prob_b1 + prob_d1 + prob_b2) && r > (prob_b1 + prob_d1)) {

                                    n_b_seg2++;
                                } else if(r > (prob_b1 + prob_d1 + prob_b2) && r <= (prob_b1 + prob_d1 + prob_b2 + prob_d2)) {
                                    n_d_seg2++;
                                }

                            }
//                                Integer b_mut_seg1 = Poisson.staticNextInt(params.rawMutationRate * 0.5 * params.p_b_seg1); //double check if rate vs mean
//                                Integer d_mut_seg1 = (int)Math.round(mut1 * params.p_d_seg1);
//                                Integer b_mut_seg2 = (int)Math.round(mut2 * params.p_b_seg2);
//                                Integer d_mut_seg2 = (int)Math.round(mut2 * params.p_d_seg2);


                            //System.out.println((mut1+mut2)+","+(b_mut_seg1+b_mut_seg2+d_mut_seg1+d_mut_seg2));

                            double seg1fitness_new = 0.0;
                            double seg2fitness_new = 0.0;


                            for (int s1 = 0; s1 < n_b_seg1; s1++) {

                                double sb = params.benDFE.nextDouble();

                                if (sb < 0 || sb > 1) {

                                    while (sb < 0 || sb > 1) {

                                        sb = params.benDFE.nextDouble();
                                        //System.out.println("sb "+sb);

                                    }
                                }


                                seg1fitness_new += (Math.log10(1 + sb));

                            }

                            for (int s1 = 0; s1 < n_d_seg1; s1++) {

                                double sd = params.delDFE.nextDouble();

                                if (sd < 0 || sd > 1) {

                                    while (sd < 0 || sd > 1) {

                                        sd = params.delDFE.nextDouble();

                                    }
                                }

                                seg1fitness_new += (Math.log10(1 - sd));
                            }

                            for (int s2 = 0; s2 < n_b_seg2; s2++) {

                                double sb = params.benDFE.nextDouble();

                                if (sb < 0 || sb > 1) {

                                    while (sb < 0 || sb > 1) {

                                        sb = params.benDFE.nextDouble();


                                    }
                                }

                                seg2fitness_new += (Math.log10(1 + sb));
                            }

                            for (int s2 = 0; s2 < n_d_seg2; s2++) {

                                double sd = params.delDFE.nextDouble();//params.delDFE.nextDouble();

                                if (sd < 0 || sd > 1) {

                                    while (sd < 0 || sd > 1) {

                                        sd = params.delDFE.nextDouble();//params.delDFE.nextDouble();

                                    }
                                }

                                seg2fitness_new += (Math.log10(1 - sd));
                            }


                            newFitness += (seg1fitness_new + seg2fitness_new);


                            iMatrix.logId(I_matrix_length);
                            iMatrix.logBirth(t_curr);
                            iMatrix.logDeath(Double.NEGATIVE_INFINITY);
                            iMatrix.logParent(parent);
                            iMatrix.logParentCo(parent_co);
                            iMatrix.logReassortant(reassortant);
                            iMatrix.logFitness(newFitness);
                            iMatrix.logSegFitness(seg1Fitness + seg1fitness_new);
                            iMatrix.logSeg1parent(seg1parent);
                            iMatrix.logSeg2parent(seg2parent);


                            I_matrix_curr.add(I_matrix_length);
                            I_matrix_length++;


                            j++;


                        }
                        break;

                }
                event++;
            }
            t_curr = t_next;


        }
        System.out.println(">" + I_matrix_length);

        double meanTMRCA1 = stats1.getMean();
        double meanTMRCA2 = stats2.getMean();
        double meanDiff = stats3.getMean();

        double std = stats1.getStandardDeviation();
        double min = stats1.getMin();
        double max = stats1.getMax();
        double lowerq = stats1.getPercentile(25);
        double upperq = stats1.getPercentile(75);
        double median1 = stats1.getPercentile(50);
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

        System.out.println("median1: " + median1);
        System.out.println("std: " + std);
        System.out.println("min: " + min + ", max: " + max);
        System.out.println("LQ: " + lowerq + ", UQ: " + upperq);
        System.out.println();
        System.out.println("median2: " + median2);
        System.out.println("std: " + std2);
        System.out.println("min: " + min2 + ", max: " + max2);
        System.out.println("LQ: " + lowerq2 + ", UQ: " + upperq2);
        System.out.println();
        System.out.println("median3: " + median3);
        System.out.println("std: " + std3);
        System.out.println("min: " + min3 + ", max: " + max3);
        System.out.println("LQ: " + lowerq3 + ", UQ: " + upperq3);

        System.out.println();
        System.out.println("med_var: "+var.getPercentile(50)+", geo_mean: "+var.getGeometricMean());


        Set<Double> generations = fitDistTime.keySet();

        List<Double> sorted_generations = new ArrayList<Double>();

        sorted_generations.addAll(generations);

        Collections.sort(sorted_generations);




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

    public double getSegPopFitness(List<Integer> I_curr, int seg) {

        Double totalFitness = 0.0;
        int count = 0;
        for(Integer i: I_curr) {

            if(seg==0) {
                if(iMatrix.getSegFitness(i) > 0) {
                    totalFitness += iMatrix.getSegFitness(i);
                    count++;
                }
            }
            else {
                if (iMatrix.getFitness(i) > 0 && iMatrix.getSegFitness(i) > 0) {

                    totalFitness += iMatrix.getFitness(i)-iMatrix.getSegFitness(i);
                    count++;
                }
            }
        }


        return (totalFitness/count);


    }


    public double getPopFitnessVar(double popFitness, List<Integer> I_curr) {

        Double totalSquaredDeviations = 0.0;

        int count = 0;
        for(Integer i: I_curr) {

            double fitness = iMatrix.getFitness(i);
            if(fitness != Double.NEGATIVE_INFINITY) {
                totalSquaredDeviations += Math.pow((fitness - popFitness), 2);
            }
            else{
                count++;
            }
        }

        return (totalSquaredDeviations)/(I_curr.size()-count);

    }

    public double getPopFitnessVar(double segFitness, List<Integer> I_curr, int seg) {

        Double totalSquaredDeviations = 0.0;

        int count = 0;

        for(Integer i: I_curr) {

            if(seg==0) {
                double fitness = iMatrix.getSegFitness(i);
                if (fitness != Double.NEGATIVE_INFINITY) {
                    totalSquaredDeviations += Math.pow((fitness - segFitness), 2);
                } else {
                    count++;
                }
            }
            else{
                double fitness = iMatrix.getFitness(i)-iMatrix.getSegFitness(i);
                if (fitness != Double.NEGATIVE_INFINITY) {
                    totalSquaredDeviations += Math.pow((fitness - segFitness), 2);
                } else {
                    count++;
                }

            }
        }

        return (totalSquaredDeviations)/(I_curr.size()-count);

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



        List<Integer> indices1 = new ArrayList<Integer>();
        List<Integer> indices2 = new ArrayList<Integer>();

        for(int i= 0; i < I_matrix_curr.size(); i++) {

            indices1.add(i);
            indices2.add(i);
        }

        Collections.shuffle(indices1, new Random((long)(Math.floor(Math.random()*10000))));
        Collections.shuffle(indices2, new Random((long)(Math.floor(Math.random()*10000))));
        int samplingDepth = I_matrix_curr.size()/100;

        for (int i = 0; i < 100; i++) {


            int j = samplingDepth*i;
            //System.out.println(j);
            Integer vA = I_matrix_curr.get(i);
            //int r = (int)Math.floor(Math.random()*I_matrix_curr.size());
            Integer vB = I_matrix_curr.get(indices2.get(i));


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


    @Override
    public double getRatesSum(DoubleArrayList rates) {
        return 0;
    }

    @Override
    public void runSimulation(EpiParams params) {

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

    @Override
    public List<Integer> getCurrentInfected(int Y_curr, int startingNumber) {

        List<Integer> I_matrix_curr = new ArrayList<Integer>();
        int i = startingNumber;
        while(i<Y_curr){

            I_matrix_curr.add(i);
            i++;
        }

        return I_matrix_curr;

    }

    @Override
    public infectionHistory getIMatrix() {
        return iMatrix;
    }


    @Override
    public coinfectionHistory getIcoMatrix() {
        return icoMatrix;
    }

    private Comparator<String> stringComparator = new Comparator<String>() {
        @Override
        public int compare(String o1, String o2) {
            double diff = Double.parseDouble(o1)-Double.parseDouble(o2);
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

    private  Comparator<Integer>  comparator1 = new Comparator<Integer>() {
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

    private  Comparator<Integer>  comparator2 = new Comparator<Integer>() {
        @Override
        public int compare(Integer o1, Integer o2) {




            int icoParentA1 = icoMatrix.getParent1(o1);
            int icoParentA2 = icoMatrix.getParent2(o1);

            double meanIcoA_fitness = (iMatrix.getFitness(icoParentA1)+iMatrix.getFitness(icoParentA2))/2;

            int icoParentB1 = icoMatrix.getParent1(o2);
            int icoParentB2 = icoMatrix.getParent2(o2);

            double meanIcoB_fitness = (iMatrix.getFitness(icoParentB1)+iMatrix.getFitness(icoParentB2))/2;

            double diff = meanIcoA_fitness - meanIcoB_fitness;

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
    private List<Double> calculateWeights(List<Integer> I_matrix_curr) {

        List<Double> cumulativeWeights = new ArrayList<Double>();
        List<Double> fitnessValues = new ArrayList<Double>();
        List<Double> relativeFitnessValues = new ArrayList<Double>();

        Collections.sort(I_matrix_curr, comparator1);


        for(Integer i: I_matrix_curr) {

            fitnessValues.add(iMatrix.getFitness(i));
        }

        //Double minFit = Collections.min(fitnessValues);
        Double maxFit = Collections.max(fitnessValues);


        for(Integer i: I_matrix_curr) {

            Double fitness = iMatrix.getFitness(i);


            Double relfitness = fitness/maxFit;//(Math.pow(10, fitness-maxFit));//Math.log10(Math.pow(10, fitness/maxFit));//(fitness/maxFit)-1.0);


            relativeFitnessValues.add(relfitness);
        }
        //int count = Collections.frequency(relativeFitnessValues,0.0);

        return relativeFitnessValues;
        //return relativeFitnessValues;

    }

    private List<Double> calculateIcoWeights(List<Integer> Ico_matrix_curr) {

        List<Double> cumulativeWeights = new ArrayList<Double>();
        List<Double> fitnessValues = new ArrayList<Double>();
        List<Double> relativeFitnessValues = new ArrayList<Double>();

        Collections.sort(Ico_matrix_curr, comparator2);


        for(Integer i: Ico_matrix_curr) {

            int icoParent1 = icoMatrix.getParent1(i);
            int icoParent2 = icoMatrix.getParent2(i);

            double meanIco_fitness = (iMatrix.getFitness(icoParent1)+iMatrix.getFitness(icoParent2))/2;


            fitnessValues.add(meanIco_fitness);
        }

        //Double minFit = Collections.min(fitnessValues);
        Double maxFit = Collections.max(fitnessValues);


        for(Integer i: Ico_matrix_curr) {


            int icoParent1 = icoMatrix.getParent1(i);
            int icoParent2 = icoMatrix.getParent2(i);

            double meanIco_fitness = (iMatrix.getFitness(icoParent1)+iMatrix.getFitness(icoParent2))/2;


            //Double fitness = iMatrix.getFitness(i);


            Double relfitness = meanIco_fitness/maxFit;//(Math.pow(10, fitness-maxFit));//Math.log10(Math.pow(10, fitness/maxFit));//(fitness/maxFit)-1.0);


            relativeFitnessValues.add(relfitness);
        }
        //int count = Collections.frequency(relativeFitnessValues,0.0);

        return relativeFitnessValues;
        //return relativeFitnessValues;

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

            //Double real_fitness = iMatrix.getFitness(parent_index);
            return parent_index;

        }

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

    private void chooseDeath(double t_curr, List<Integer> I_matrix_curr, List<Integer> Ico_matrix_curr, int Yr_curr, int Yr_primary_curr, boolean Is_or_Ico) {

        if(Is_or_Ico) {

            // Is recovers

            int i2 = (int) Math.floor(Math.random() * I_matrix_curr.size());

            Integer dead_Is = I_matrix_curr.get(i2);
            I_matrix_curr.remove(i2);

            int reassortant = iMatrix.getReassortant(dead_Is);
            int parent_co = iMatrix.getParentCo(dead_Is);

            if (reassortant == 1) {

                Yr_curr--;
                //cumulative Irs_i (secondary Ir)

                if (parent_co >= 0 && Yr_primary_curr > 0) {
                    //System.out.println(">"+Yr_primary_curr);
                    Yr_primary_curr--;
                    //System.out.println(Yr_primary_curr);

                }
            }

            iMatrix.setDeath(dead_Is, t_curr);
        }
        else{

            //Ico recovers
            int i2 = (int) Math.floor(Math.random() * Ico_matrix_curr.size());

            Integer dead_Ico = Ico_matrix_curr.get(i2);
            Ico_matrix_curr.remove(i2);

            icoMatrix.setDeath(dead_Ico, t_curr);
        }

    }

    public static void main(String [] args) {

        EpiParams params = new EpiParams();

        params.p_b_seg1 = Double.parseDouble(args[0]);
        params.p_b_seg2 = Double.parseDouble(args[1]);
        params.p_d_seg1 = Double.parseDouble(args[2]);
        params.p_d_seg2 = Double.parseDouble(args[3]);
        params.rawMutationRate = Double.parseDouble(args[4]);
        params.psi = Double.parseDouble(args[5]);
        int segs = Integer.parseInt(args[6]);
        params.relativeRate = Double.parseDouble(args[7]);
        params.Is_init = Integer.parseInt(args[8]);
        params.Ico_init = Integer.parseInt(args[9]);
        params.s_ben = Double.parseDouble(args[10]);
        params.s_del = Double.parseDouble(args[11]);

        int no_of_sims = Integer.parseInt(args[12]);
        boolean getTrees = Boolean.parseBoolean(args[13]);

        double[][] tmrca_seg1 = new double[no_of_sims][61];
        double[][] tmrca_seg2 = new double[no_of_sims][61];
        double[][] fitness = new double[no_of_sims][61];
        double[][] fitnessv = new double[no_of_sims][61];
        double[][] fitness_1 = new double[no_of_sims][61];
        double[][] fitnessv_1 = new double[no_of_sims][61];
        double[][] fitness_2 = new double[no_of_sims][61];
        double[][] fitnessv_2 = new double[no_of_sims][61];

        params.benDFE = new Exponential(1/params.s_ben, params.randomGenerator);
        params.delDFE = new Exponential(1/params.s_del, params.randomGenerator);



        params.print();

        for(int i=0; i < no_of_sims; i++) {
            System.out.println("sim "+(i+1));
            moranReassortment model = new moranReassortment();
            //reassortmentTauLeap model = new reassortmentTauLeap();


            model.runSimulation(params, i, tmrca_seg1, tmrca_seg2, fitness, fitnessv, fitness_1, fitnessv_1, fitness_2, fitnessv_2);

            SimulateTree tree = new SimulateTree();

//            tree.n_lineages = 100;
//
            if(getTrees) {
                tree.sampleStartTime = 40 * 365;
                tree.sampleEndTime = params.simulationEndTime;
                tree.n_lineages = 300;
                //System.out.println(tree.samplingSchemeForOnePatchModel);

                String treeFilename = "tree_N_"+tree.n_lineages+"N_"+ params.rawMutationRate+"_sb_"+params.s_ben+"_sd_"+params.s_del+"_pb1_"+params.p_b_seg1+"_pb2_"+params.p_b_seg2+"_pd1_"+params.p_d_seg1+"_pd2_"+params.p_d_seg2+"_psi_"+params.psi+"_co_"+params.relativeRate+"_N_"+(params.Is_init+params.Ico_init)+"_simTime_"+params.simulationTime+"yrs_simNo_"+(i+1)+"_simNo_"+(i+1)+"_segment_";

                tree.getTransmissionTrees(segs, model, treeFilename, (i+1));
            }


        }



        File tmrca_seg1_output = new File("tmrca1_through_time_mu_"+params.rawMutationRate+"_sb_"+params.s_ben+"_sd_"+params.s_del+"_pb1_"+params.p_b_seg1+"_pb2_"+params.p_b_seg2+"_pd1_"+params.p_d_seg1+"_pd2_"+params.p_d_seg2+"_psi_"+params.psi+"_co_"+params.relativeRate+"_D_"+params.durationOfInfection+"_N_"+(params.Is_init+params.Ico_init)+"_simTime_"+params.simulationTime+"yrs_nsimsz_"+no_of_sims+".csv");
        File tmrca_seg2_output = new File("tmrca2_through_time_mu_"+params.rawMutationRate+"_sb_"+params.s_ben+"_sd_"+params.s_del+"_pb1_"+params.p_b_seg1+"_pb2_"+params.p_b_seg2+"_pd1_"+params.p_d_seg1+"_pd2_"+params.p_d_seg2+"_psi_"+params.psi+"_co_"+params.relativeRate+"_D_"+params.durationOfInfection+"_N_"+(params.Is_init+params.Ico_init)+"_simTime_"+params.simulationTime+"yrs_nsimsz_"+no_of_sims+".csv");
        //File diff_in_tmrca = new File("diff_in_tmrca_through_time_mu_"+params.rawMutationRate+"_sb_"+params.s_ben+"_sd_"+params.s_del+"_pb1_"+params.p_b_seg1+"_pb2_"+params.p_b_seg2+"_pd1_"+params.p_d_seg1+"_pd2_"+params.p_d_seg2+"_psi_"+params.psi+"_co_"+params.relativeRate+"_D_"+params.durationOfInfection+"_N_"+(params.Is_init+params.Ico_init)+"_simTime_"+params.simulationTime+"yrs_nsimsz_"+no_of_sims+".csv");

        File fitness_output = new File("fitness_through_time_mu_"+params.rawMutationRate+"_sb_"+params.s_ben+"_sd_"+params.s_del+"_pb1_"+params.p_b_seg1+"_pb2_"+params.p_b_seg2+"_pd1_"+params.p_d_seg1+"_pd2_"+params.p_d_seg2+"_psi_"+params.psi+"_co_"+params.relativeRate+"_D_"+params.durationOfInfection+"_N_"+(params.Is_init+params.Ico_init)+"_simTime_"+params.simulationTime+"yrs_nsimsz_"+no_of_sims+".csv");
        File fitness_output_1 = new File("seg1_fitness_through_time_mu_"+params.rawMutationRate+"_sb_"+params.s_ben+"_sd_"+params.s_del+"_pb1_"+params.p_b_seg1+"_pb2_"+params.p_b_seg2+"_pd1_"+params.p_d_seg1+"_pd2_"+params.p_d_seg2+"_psi_"+params.psi+"_co_"+params.relativeRate+"_D_"+params.durationOfInfection+"_N_"+(params.Is_init+params.Ico_init)+"_simTime_"+params.simulationTime+"yrs_nsimsz_"+no_of_sims+".csv");
        File fitness_output_2 = new File("seg2_fitness_through_time_mu_"+params.rawMutationRate+"_sb_"+params.s_ben+"_sd_"+params.s_del+"_pb1_"+params.p_b_seg1+"_pb2_"+params.p_b_seg2+"_pd1_"+params.p_d_seg1+"_pd2_"+params.p_d_seg2+"_psi_"+params.psi+"_co_"+params.relativeRate+"_D_"+params.durationOfInfection+"_N_"+(params.Is_init+params.Ico_init)+"_simTime_"+params.simulationTime+"yrs_nsimsz_"+no_of_sims+".csv");


        try {
            BufferedWriter writer1 = new BufferedWriter(new FileWriter(tmrca_seg1_output));
            BufferedWriter writer2 = new BufferedWriter(new FileWriter(tmrca_seg2_output));
            BufferedWriter writer3 = new BufferedWriter(new FileWriter(fitness_output));
            //BufferedWriter writer4 = new BufferedWriter(new FileWriter(diff_in_tmrca));
            BufferedWriter writer5 = new BufferedWriter(new FileWriter(fitness_output_1));
            BufferedWriter writer6 = new BufferedWriter(new FileWriter(fitness_output_2));

            for(int i=0; i < 61; i++) {

                String tmrca1_time = "";
                String tmrca1_raw = "";
                String tmrca2_time = "";
                String tmrca2_raw = "";
                String fitness_string = "";
                String fitnessv_string = "";
                String tmrca_diff = "";
                String fitness_string_1 = "";
                String fitnessv_string_1 = "";
                String fitness_string_2 = "";
                String fitnessv_string_2 = "";


                for (int j = 0; j < no_of_sims; j++) {

                    tmrca1_raw += tmrca_seg1[j][i] + ",";
                    tmrca2_raw += tmrca_seg2[j][i] + ",";
                    tmrca1_time += (60  - tmrca_seg1[j][i]) + ",";
                    tmrca2_time += (60 - tmrca_seg2[j][i]) + ",";
                    fitness_string += fitness[j][i] + ",";
                    fitnessv_string += fitnessv[j][i] + ",";
                    //tmrca_diff += tmrca_seg1[j][i]-tmrca_seg2[j][i] + ",";
                    fitness_string_1 += fitness_1[j][i] + ",";
                    fitnessv_string_1 += fitnessv_1[j][i] + ",";
                    fitness_string_2 += fitness_2[j][i] + ",";
                    fitnessv_string_2 += fitnessv_2[j][i] + ",";

                }

                writer1.write(i + "," + tmrca1_time+tmrca1_raw + "\n");
                writer2.write(i + "," + tmrca2_time+tmrca2_raw + "\n");
                writer3.write(i + "," + fitness_string+fitnessv_string + "\n");
                writer5.write(i + "," + fitness_string_1+fitnessv_string_1 + "\n");
                writer6.write(i + "," + fitness_string_2+fitnessv_string_2 + "\n");

                //writer4.write(i + "," + tmrca_diff + "\n");

                writer1.flush();
                writer2.flush();
                writer3.flush();
               // writer4.flush();
                writer5.flush();
                writer6.flush();

            }

            writer1.close();
            writer2.close();
            writer3.close();
            //writer4.close();
            writer5.close();
            writer6.close();


        } catch (IOException e) {
            e.printStackTrace();
        }



    }


}
