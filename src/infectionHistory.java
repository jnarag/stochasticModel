import java.util.ArrayList;
import java.util.List;

/**
 * Created with IntelliJ IDEA.
 * User: Jayna
 * Date: 25/11/2012
 * Time: 11:02
 * To change this template use File | Settings | File Templates.
 */
public class infectionHistory implements iHistory {

    List<Integer> id = null;
    List<Double> birth = null;
    List<Double> death = null;
    List<Integer> parent = null;
    List<Integer> parentCo = null;
    List<Integer> reassortant = null;
    List<Double> fitness = null;
    List<Double> segment1Fitness = null;
    List<Integer> patch = null;

    public infectionHistory (int size) {

        id = new ArrayList<Integer>(size);
        birth = new ArrayList<Double>(size);
        death = new ArrayList<Double>(size);
        parent = new ArrayList<Integer>(size);
        parentCo = new ArrayList<Integer>(size);
        reassortant = new ArrayList<Integer>(size);
        fitness = new ArrayList<Double>(size);
        segment1Fitness = new ArrayList<Double>(size);
        patch = new ArrayList<Integer>(size);

    }

    public infectionHistory () {

        id = new ArrayList<Integer>();
        birth = new ArrayList<Double>();
        death = new ArrayList<Double>();
        parent = new ArrayList<Integer>();
        parentCo = new ArrayList<Integer>();
        reassortant = new ArrayList<Integer>();
        fitness = new ArrayList<Double>();
        segment1Fitness = new ArrayList<Double>();
        patch = new ArrayList<Integer>();

    }

    public infectionHistory (infectionHistory iMatrix) {

        id = new ArrayList<Integer>();
        birth = new ArrayList<Double>();
        death = new ArrayList<Double>();
        parent = new ArrayList<Integer>();
        parentCo = new ArrayList<Integer>();
        reassortant = new ArrayList<Integer>();
        fitness = new ArrayList<Double>();
        segment1Fitness = new ArrayList<Double>();
        patch = new ArrayList<Integer>();

        id.addAll(iMatrix.getId());
        birth.addAll(iMatrix.getBirth());
        death.addAll(iMatrix.getDeath());
        parent.addAll(iMatrix.getParent());
        parentCo.addAll(iMatrix.getParentCo());
        reassortant.addAll(iMatrix.getReassortant());
        fitness.addAll(iMatrix.getFitness());
        segment1Fitness.addAll(iMatrix.getSegFitness());
        patch.addAll(iMatrix.getPatch());

    }

    public infectionHistory (List<Integer> id, List<Double> birth, List<Double> death) {

        this.id = id;
        this.birth = birth;
        this.death = death;
    }

    public void logId(int id) {

        this.id.add(id);
    }
    public void logBirth(double birth){

        this.birth.add(birth);
    }
    public void logDeath(double death) {

        this.death.add(death);

    }
    public void logParent(int parent) {

        this.parent.add(parent);
    }
    public void logParentCo(int parentCo) {

        this.parentCo.add(parentCo);
    }
    public void logFitness(double fitness) {

        this.fitness.add(fitness);
    }
    public void logSegFitness(double fitness) {

        this.segment1Fitness.add(fitness);
    }
    public void logReassortant(int reassortant) {

        this.reassortant.add(reassortant);
    }
    public void logPatch(int patch) {

        this.patch.add(patch);
    }

    public void setBirth(int index, double birth){

        this.birth.set(index, birth);
    }
    public void setDeath(int index, double death) {

        this.death.set(index, death);

    }
    public void setParent(int index, Integer parent) {

        this.parent.set(index, parent);
    }
    public void setParentCo(int index, Integer parentCo) {

        this.parentCo.set(index, parentCo);
    }
    public void setFitness(int index, double fitness) {

        this.fitness.set(index, fitness);
    }
    public void setSegFitness(int index, double fitness) {

        this.segment1Fitness.set(index, fitness);
    }
    public void setReassortant(int index, int reassortant) {

        this.reassortant.set(index, reassortant);
    }
    public void setPatch(int index, int patch) {

        this.patch.set(index, patch);
    }

    public int getId(int index) {

        return this.id.get(index);
    }
    public double getBirth(int index){

        return this.birth.get(index);
    }
    public double getDeath(int index) {

        return this.death.get(index);

    }
    public int getParent(int index) {

        return this.parent.get(index);
    }
    public int getParentCo(int index) {

        return this.parentCo.get(index);
    }
    public double getFitness(int index) {

        return this.fitness.get(index);

    }
    public double getSegFitness(int index) {

        return this.segment1Fitness.get(index);

    }
    public int getReassortant(int index) {

        return this.reassortant.get(index);
    }
    public int getPatch(int index) {

        return this.patch.get(index);
    }

    public void addEntry(infectionHistory matrix, int i) {

        logId(matrix.getId(i));
        logBirth(matrix.getBirth(i));
        logDeath(matrix.getDeath(i));
        logParent(matrix.getParent(i));
        logParentCo(matrix.getParentCo(i));
        logReassortant(matrix.getReassortant(i));
        logFitness(matrix.getFitness(i));
        logSegFitness(matrix.getSegFitness(i));
        logPatch(matrix.getPatch(i));
    }

    public List<Integer> getId() {
        return this.id;
    }
    public List<Double> getBirth() {
        return this.birth;
    }
    public List<Double> getDeath() {
        return this.death;
    }
    public List<Integer> getParent() {
        return this.parent;
    }
    public List<Integer> getParentCo() {
        return this.parentCo;
    }
    public List<Integer> getReassortant() {
        return this.reassortant;
    }
    public List<Double> getFitness() {
        return this.fitness;
    }
    public List<Double> getSegFitness() {
        return this.segment1Fitness;
    }
    public List<Integer> getPatch() {
        return this.patch;
    }


    public infectionHistory getInfectionHistory() {

       infectionHistory iMatrix = new infectionHistory();


        iMatrix.id.addAll(this.id);
        iMatrix.birth.addAll(this.birth);
        iMatrix.death.addAll(this.death);
        iMatrix.parent.addAll(this.parent);
        iMatrix.parentCo.addAll(this.parentCo);
        iMatrix.reassortant.addAll(this.reassortant);
        iMatrix.fitness.addAll(this.fitness);
        iMatrix.segment1Fitness.addAll(this.segment1Fitness);
        iMatrix.patch.addAll(this.patch);

        return iMatrix;
    }

}
