/**
 * Created with IntelliJ IDEA.
 * User: jr207
 * Date: 2/7/13
 * Time: 4:02 PM
 * To change this template use File | Settings | File Templates.
 */
public interface iHistory {

    void logBirth(double birth);
    void logDeath(double death);
    //void logParent(int parent);
    //void logParentCo(int parentCo);
    void logFitness(double fitness);
    void logPatch(int patch);
    void setBirth(int index, double birth);
    void setDeath(int index, double death);
    //void setParent(int index, Integer parent);
    //void setParentCo(int index, Integer parent);
    void setFitness(int index, double fitness);
    //void setReassortant(int index, int reassortant);
    void setPatch(int index, int reassortant);
    double getBirth(int index);
    double getDeath(int index);
    //int getParent(int index);
    //void getParentCo(int index);
    double getFitness(int index);
    //void getReassortant(int index);
    int getPatch(int index);


}
