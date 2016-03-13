import java.util.ArrayList;
import java.util.List;

/**
 * Created with IntelliJ IDEA.
 * User: Jayna
 * Date: 25/11/2012
 * Time: 11:02
 * To change this template use File | Settings | File Templates.
 */
public class coinfectionHistory implements iHistory {


    List<Double> birth = null;
    List<Double> death = null;
    List<Integer> parent1 = null;
    List<Integer> parent2 = null;
    List<Double> fitness = null;
    List<Integer> patch = null;

    public void initialize(int size) {

        birth = new ArrayList<Double>(size);
        death = new ArrayList<Double>(size);
        parent1 = new ArrayList<Integer>(size);
        parent2 = new ArrayList<Integer>(size);
        //fitness = new ArrayList<Double>(size);
        patch = new ArrayList<Integer>(size);

    }

    public void logBirth(double birth){

        this.birth.add(birth);
    }

    public void logDeath(double death) {

        this.death.add(death);

    }

    public void logParent1(int parent1) {

        this.parent1.add(parent1);
    }

    public void logParent2(int parent2) {

        this.parent2.add(parent2);
    }

    public void logFitness(double fitness) {

        this.fitness.add(fitness);
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

    public void setParent1(int index, int parent1) {

        this.parent1.set(index, parent1);
    }

    public void setParent2(int index, int parent2) {

        this.parent2.set(index, parent2);
    }

    public void setFitness(int index, double fitness) {

        this.fitness.set(index, fitness);
    }

    public void setPatch(int index, int patch) {

        this.patch.set(index, patch);
    }

    public double getBirth(int index){

        return this.birth.get(index);
    }

    public double getDeath(int index) {

        return this.death.get(index);

    }

    public int getParent1(int index) {

        return this.parent1.get(index);
    }

    public int getParent2(int index) {

        return this.parent2.get(index);
    }

    public double getFitness(int index) {

        return this.fitness.get(index);
    }

    public int getPatch(int index) {

        return this.patch.get(index);
    }

    public coinfectionHistory getCoinfectionHistory(){


        coinfectionHistory icoMatrix = new coinfectionHistory();


        icoMatrix.birth.addAll(this.birth);
        icoMatrix.death.addAll(this.death);
        icoMatrix.parent1.addAll(this.parent1);
        icoMatrix.parent2.addAll(this.parent2);
        icoMatrix.fitness.addAll(this.fitness);
        icoMatrix.patch.addAll(this.patch);

        return icoMatrix;

    }
}
