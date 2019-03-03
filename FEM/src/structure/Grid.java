package structure;
import matrixEquation.MatrixEquation;

import java.math.RoundingMode;
import java.text.DecimalFormat;
import java.util.HashMap;
import java.util.Map;


public class Grid {
    private Node [][] nodes;
    private Element[][] elements;
    private Map<Integer,Node> nodesMap=new HashMap<>();
    private double[][] hMatrix;
    private double[][] cMatrix;
    private double [] pVector;
    private double [] currentTemperature;
    private double H;
    private double L;
    private int nH;
    private int nL;

    public Grid(double H, double L, int nH, int nL ) {
        nodes=new Node[nL][nH];
        elements=new Element[nL-1][nH-1];
        hMatrix=new double [nH*nL][nH*nL];
        cMatrix=new double [nH*nL][nH*nL];
        pVector=new double [nH*nL];
        currentTemperature =new double [nH*nL];

        this.H=H;
        this.L=L;
        this.nH=nH;
        this.nL=nL;

        //INITIATE NODES
        double lStep=L/(nL-1);
        double hStep=H/(nH-1);

        for(int i=0;i<nL;++i) {
            for(int j=0;j<nH;++j) {
                nodes[i][j]=new Node();
                nodes[i][j].setId(i*nH+j+1);
                nodes[i][j].setX(i*lStep);
                nodes[i][j].setY(j*hStep);
                //nodes[i][j].setT0(t0);//////////////
                nodesMap.put((i*nH+j+1),nodes[i][j]);
            }
        }
        //INITIATE ELEMENTS
        for(int i=0;i<(nL-1);++i){
            for (int j=0;j<(nH-1);++j){
                elements[i][j]=new Element();
                elements[i][j].setId(i*nH+j+1);
                elements[i][j].setNodes(i*nH+j+1,i*nH+j+1+nH,i*nH+j+1+nH+1,i*nH+j+1+1);
            }
        }

    }

    public void countGrid(double k, double c, double ro, double alfa,double tEnv, double deltaTime, double [] initialTemperature  ){
        boolean[]bc=new boolean[4];//boundary conditions
        double tmpMatrix[][]=new double[nH*nL][nH*nL];
        double tmpVector[]=new double [nH*nL];

        for(int i=0;i<(nL-1);++i) {//for all elements in the grid:
            for (int j = 0; j < (nH - 1); ++j) {

                if(j==0){ bc[0]=true;}
                else {bc[0]=false; }
                if(i==nL-2){bc[1]=true;}
                else {bc[1]=false; }
                if(j==nH-2){bc[2]=true;}
                else {bc[2]=false; }
                if(i==0){bc[3]=true;}
                else {bc[3]=false; }

                elements[i][j].countElement(nodesMap.get(i*nH+j+1),nodesMap.get(i*nH+j+1+nH),nodesMap.get(i*nH+j+1+nH+1),nodesMap.get(i*nH+j+1+1),bc, k,c,ro,alfa,tEnv);
                for(int a=0;a<4;++a){//creating global matrix from local
                    for (int b=0;b<4;++b){
                        hMatrix[(elements[i][j].nodes[a])-1][(elements[i][j].nodes[b])-1]+=elements[i][j].hMatrix[a][b];
                        cMatrix[(elements[i][j].nodes[a])-1][(elements[i][j].nodes[b])-1]+=elements[i][j].cMatrix[a][b];
                    }
                    pVector[(elements[i][j].nodes[a])-1]+=elements[i][j].pVector[a];
                }
            }
        }

        for(int i=0;i<(nH*nL);++i) {//for all elements in the h/cMatrix:!!!!!!!!!!!!!!!!!!!!!!!!!!!
            for (int j = 0; j < (nH*nL); ++j) {
                tmpMatrix[i][j]=hMatrix[i][j]+cMatrix[i][j]/deltaTime;
                tmpVector[i]+=(cMatrix[i][j]/deltaTime)*initialTemperature[j];
            }
            tmpVector[i]+=pVector[i];
        }
        MatrixEquation equation=new MatrixEquation();

        currentTemperature= equation.solveEquationGauss(tmpMatrix,tmpVector);
//        System.out.println("Current temperature vector: ");
//        for(int j=0; j<nL*nH;++j){
//            System.out.println(currentTemperature[j]+" ");
//        }
//        System.out.println();

    }

    public void showGlobalH(){

        DecimalFormat df = new DecimalFormat("#.####");
        df.setRoundingMode(RoundingMode.CEILING);

        System.out.println("Global H matrix: ");
        for(int i=0; i<nL*nH;++i){
            for(int j=0; j<nL*nH;++j){
                System.out.print(df.format(hMatrix[i][j])+" ");
            }
            System.out.println();
        }
        System.out.println();
    }
    public void showGlobalC(){
        DecimalFormat df = new DecimalFormat("#.####");
        df.setRoundingMode(RoundingMode.CEILING);

        System.out.println("Global C matrix: ");
        for(int i=0; i<nL*nH;++i){
            for(int j=0; j<nL*nH;++j){
                System.out.print(df.format(cMatrix[i][j])+" ");
            }
            System.out.println();
        }
        System.out.println();
    }
    public void showGlobalP(){
        DecimalFormat df = new DecimalFormat("#.####");
        df.setRoundingMode(RoundingMode.CEILING);

        System.out.println("Global P vector: ");
        for(int j=0; j<nL*nH;++j){
            System.out.print(df.format(pVector[j])+" ");
        }
        System.out.println();
    }
    public void showTemperature(){

        DecimalFormat df = new DecimalFormat("#.####");
        df.setRoundingMode(RoundingMode.CEILING);

        System.out.println("Current temperature vector: ");
        for(int j=0; j<nL*nH;++j){
            System.out.print(df.format(currentTemperature[j])+" ");
        }
        System.out.println();
    }

    public double getMinTemperature(){
        double min=currentTemperature[0];
        for(int j=1; j<nL*nH;++j){
           if(currentTemperature[j]<min){
               min=currentTemperature[j];
           }
        }
        return min;
    }
    public double getMaxTemperature(){
        double max=currentTemperature[0];
        for(int j=1; j<nL*nH;++j){
            if(currentTemperature[j]>max){
                max=currentTemperature[j];
            }
        }
        return max;
    }
    public double[] getCurrentTemperature() {
        return currentTemperature;
    }
}
