import cern.colt.list.DoubleArrayList;

import java.util.List;

/**
 * Created with IntelliJ IDEA.
 * User: jr207
 * Date: 8/28/12
 * Time: 5:24 PM
 * To change this template use File | Settings | File Templates.
 */
public interface EpiModel {

     infectionHistory iMatrix = new infectionHistory((int)20e6);
     coinfectionHistory icoMatrix = new coinfectionHistory();

    void runSimulation(EpiParams params);

    double getRatesSum(DoubleArrayList rates);

    int chooseEvent(DoubleArrayList rates, double randomNo, double rates_sum);

    double[] normCumSum(DoubleArrayList rates, double rateMatrix_sum);

    List<Integer> getCurrentInfected(int I_curr, int startingNumber);

    infectionHistory getIMatrix();

    coinfectionHistory getIcoMatrix();

    //void getCoinfectedHistory(int[] coinfectHistory, int birth, int death, int parent1, int parent2);

}
