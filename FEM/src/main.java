import matrixEquation.MatrixEquation;
import structure.Grid;

import java.math.RoundingMode;
import java.text.DecimalFormat;

public class main {

    public static void main(String[] args) {
        double initialTemperature = 200;
        double simulationTime=500;
        double deltaTime=50;
        double ambientTemperature = 800;
        int alfa = 300;
        double H=0.200;
        double L=0.200;
        int nH=10;
        int nL=10;
        double specificHeat=910;//c
        double conductivity=273;//k
        double density=2700;//ro
        double[] t0=new double [nH*nL];
        for(int i=0;i<nH*nL;++i){
            t0[i]=initialTemperature;
        }
        DecimalFormat df = new DecimalFormat("#.####");
        df.setRoundingMode(RoundingMode.CEILING);

        Grid grid=new Grid(H,L,nH,nL);

//        grid.countGrid(conductivity,specificHeat,density,alfa,ambientTemperature,deltaTime,t0);
//        grid.showTemperature();
        for(double i=deltaTime;i<=simulationTime;i+=deltaTime){
            grid.countGrid(conductivity,specificHeat,density,alfa,ambientTemperature,deltaTime,t0);
            System.out.println("For simulation time = "+ i + " seconds" );
            //grid.showTemperature();
            System.out.println("Minimum temperature is "+df.format(grid.getMinTemperature())+" degrees");
            System.out.println("Maximum temperature is "+df.format(grid.getMaxTemperature())+" degrees");
            System.out.println();
            t0=grid.getCurrentTemperature();
        }
//        double test= grid.elements[0][0].SHAPE_FUNCTION[0][0];
//        System.out.println(grid.nodes[0][0].getId());
//        System.out.println(grid.elements[0][0].KSI[0]);
//        grid.elements[0][0].showLocalH();
//        grid.elements[0][0].showLocalC();
//        grid.elements[0][0].showLocalP();

//        double [][] testA = new double[][] {
//                new double[] { 2,3,1 },
//                new double [] {3,-1,2 },
//                new double [] { 1, 4,-1},
//        };
//        double []testB= {1,1,2};
//         MatrixEquation equation = new MatrixEquation(testA,testB);
//         double [] result;
//         result=equation.solveEquationGauss(testA,testB);
//        System.out.println(result[0]+ " "+result[1]+" "+result[2]);
        return;
    }
}
