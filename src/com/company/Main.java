package com.company;
import java.io.BufferedWriter;
import java.io.BufferedReader;
import java.io.InputStreamReader;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.math.*;
import java.util.ArrayList;

public class Main {

    public static void main(String[] args) {

        System.out.println("\nWelcome to the project simulator.\n\nPlease wait while parameters are being loaded.\n");
        // write your code here
        int nsales = 8;
        int nproductions = 3;
        int nprods = 5;
        int nmkts = 2;
        int ntrucks = 4;
        int ntimes = 6;
        int nrawmaterials = 4;
        int nvendors = 3;
        String inputpath = File.separator+"Input"+File.separator;

        //parameters
        Object[] param = fillParameters(nproductions, nsales, ntimes, nprods, nmkts, ntrucks, nrawmaterials, nvendors);
        double                  hpermonth =            (double) param[0];
        double [][]              worktime =       (double [][]) param[1];
        double []           fixedpaycheck =         (double []) param[2];
        double []        variablepaycheck =         (double []) param[3];
        double []                      lb =         (double []) param[4];
        double []                      ub =         (double []) param[5];

        double []                  volume =         (double []) param[6];
        double [][]           transp_cost =       (double [][]) param[7];
        double []              reloc_cost =         (double []) param[8];
        double []              truck_cost =         (double []) param[9];
        double [][] activation_truck_cost =       (double [][]) param[10];
        double []            truck_volume =         (double []) param[11];

        double [][]     conversion_factor =       (double [][]) param[12];
        double [][]                    ar =       (double [][]) param[13];
        double [][]                    br =       (double [][]) param[14];

        double [][]     energy_conversion =       (double [][]) param[15];
        double []                 aenergy =         (double []) param[16];
        double []        energy_threshold =         (double []) param[17];

        double []              fixedsales =         (double []) param[18];
        double []        fixedproductions =         (double []) param[19];

        double [][]  inventory_cost_param =       (double [][]) param[20];

        double [][]                     M =       (double [][]) param[21];
        double [][]                     P =       (double [][]) param[22];
        double                          a =           (double ) param[23];
        double [][]                     Q =       (double [][]) param[24];
        double                          b =           (double ) param[25];

        try {
            BufferedReader reader = new BufferedReader(new InputStreamReader(System.in));
        System.out.println("\n------------------------------------------------------\n" +
                             "------------------- BUGelato Corp. -------------------\n" +
                             "-------- we know what you eat, but you don't----------\n" +
                             "------------------------------------------------------\n\n");
        System.out.print("   Team name:");
        String s = reader.readLine();
        System.out.println("\n   Hello, "+s+"!" );

        while(true) {
            System.out.println("\n   To run the simulation please enter 'r', to quit enter 'q'.");
            System.out.print("     Enter:");
            s = reader.readLine();

            if (s.equalsIgnoreCase("r")) {

                //input
                Object[] input = fillInput(nproductions, nsales, ntimes, nprods, nmkts, ntrucks, inputpath);
                double[][][][] x = (double[][][][]) input[0];
                double[][] p = (double[][]) input[1];
                double[][] wvar = (double[][]) input[2];
                double[] wfixed = (double[]) input[3];
                double[][] investment = (double[][]) input[4];
                double[][][][] T = (double[][][][]) input[5];
                double[][][] Ti = (double[][][]) input[6];

                //quantity to derive
                double[][][][] demand = new double[nmkts][nsales][ntimes][nprods];
                double[][] q = new double[nrawmaterials][nproductions];
                double[][] energy = new double[nproductions][ntimes];
                double[] Ttot = new double[ntrucks];
                double[] Treloc = new double[ntrucks];
                int[] delta = new int[nsales];
                int[] alpha = new int[nproductions];
                double[][][] yinv = new double[nprods][ntimes][nsales];
                double[][][][] v = new double[nmkts][ntimes][nsales][nprods];
                double[][] A = fillA(x, nprods, ntimes, nsales, nproductions);


                //todo energy and raw material parameters

                //computation
                //
                //todo make quantities integer!

                q = computeQuantities(x, conversion_factor, nsales, nproductions, nprods, ntimes, nrawmaterials);

                energy = computeEnergies(x, energy_conversion, nsales, nproductions, nprods, ntimes);

                alpha = computeAlpha(x, nsales, nproductions, nprods, ntimes);

                delta = computeDelta(x, nsales, nproductions, nprods, ntimes);


                //
                //demand

                for (int m = 0; m < nmkts; m++) {
                    for (int j = 0; j < nsales; j++) {
                        for (int k = 0; k < ntimes; k++) {
                            for (int n = 0; n < nprods; n++) {
                                demand[m][j][k][n] = computeDemand(a, b, M[j][n], P[m][n], Q[m][k], A[k][j], investment[m][k], p[m][n]);
                            }
                        }
                    }
                }


                Object[] wharehouses = computeWharehouses(x, demand, p, nmkts, nsales, nproductions, nprods, ntimes);
                yinv = (double[][][]) wharehouses[0];
                v = (double[][][][]) wharehouses[1];

                Object[] trucks = computeTrucks(Ti, ntrucks, ntimes, nproductions);
                Ttot = (double[]) trucks[0];
                Treloc = (double[]) trucks[1];

                //objective value

                double[] revs = computeRevenues(v, p, nmkts, nsales, ntimes, nprods);
                double revenues = revs[0];
                double revenues0 = revs[1];
                double revenues1 = revs[2];


                //todo invert raw material and energy cost functions
                double rawmaterialcost = rawMaterialCosts(q, ar, br, nvendors, nrawmaterials, nproductions);
                double energycost = energyCosts(energy, aenergy, energy_threshold, ntimes, nproductions);
                double workerscost = labourCost(wvar, wfixed, variablepaycheck, fixedpaycheck, nproductions, ntimes);
                double fixedsalescost = fixedSales(delta, fixedsales, nsales);
                double fixedproductionscost = fixedPruductions(alpha, fixedproductions, nproductions);
                double transportationcost = transportationCost(x, transp_cost, volume, nproductions, nsales, ntimes, nprods);
                double inventorycost = inventoryCost(yinv, inventory_cost_param, nsales, ntimes, nprods);
                double fixedtruckcost = fixedTruckCost(Ttot, truck_cost, ntrucks);
                double relocationcost = relocationCost(Treloc, reloc_cost, ntrucks);
                double trucksactivetioncost = activeTruckCost(T, activation_truck_cost, nproductions, nsales, ntimes, ntrucks);
                double investments = Sum(investment);

                double profit = revenues - rawmaterialcost - energycost - workerscost
                        - fixedproductionscost - fixedsalescost - fixedtruckcost
                        - transportationcost - trucksactivetioncost - inventorycost
                        - investments - relocationcost;
                String report =
                          "  FINAL REPORT                       \n" +
                          "_____________________________________________________________\n\n" +
                          "  PROFIT:                    " + profit + "\n";
                report += "_____________________________________________________________\n";
                report += "+ Revenues                   " + revenues + "\n";
                report += "- Raw material costs         " + rawmaterialcost + "\n";
                report += "- Energy costs               " + energycost + "\n";
                report += "- Labour costs               " + workerscost + "\n";
                report += "- Fixed Prod. Center costs   " + fixedproductionscost + "\n";
                report += "- Transport costs            " + transportationcost + "\n";
                report += "- Trucs Activation costs     " + trucksactivetioncost + "\n";
                report += "- Trucks purchase costs      " + fixedtruckcost + "\n";
                report += "- Trucks relocation costs    " + relocationcost + "\n";
                report += "- Inventory costs            " + inventorycost + "\n";
                report += "- Fixed Sales Center costs   " + fixedsalescost + "\n";
                report += "- Investments                " + investments + "\n";


                //feasibility
                Object[] check = new Object[2];
                boolean checkbool = true;

                //workers//production centers w.r.t. workers

                check = checkWorkersAndSaleCenters(alpha, x, wfixed, wvar, worktime, lb, ub, hpermonth,
                        nproductions, nsales, ntimes, nprods);
                checkbool = (boolean) check[0];
                if (!checkbool)
                    report = (String) check[1];


                //salescenters
                check = checkSales(delta);
                checkbool = (boolean) check[0];
                if (!checkbool)
                    report = (String) check[1];

                //trucks
                check = checkTrucks(x, T, truck_volume, volume,
                        nproductions, nsales, ntimes, nprods, ntrucks);
                if (!checkbool)
                    report = (String) check[1];
                //nonegativity
                check = checkNonNegativity(x, T, Ti, wfixed, wvar, p, investment, nproductions, nsales, ntimes, nprods, nmkts, ntrucks);
                checkbool = (boolean) check[0];
                if (!checkbool)
                    report = (String) check[1];


                //print report into a file called finalReport.txt
                File reportfile = new File("finalReport.txt");

                BufferedWriter writer = new BufferedWriter(new FileWriter(reportfile, false));
                writer.write(report);
                writer.flush();
                writer.close();
                System.out.println("     REPORT: "+report);
            }else if(s.equalsIgnoreCase("q")){
                System.out.println("   Thank you for using the simulation tool.\n   See you next time and remember:" +
                        (Math.random()<0.3? " better some bugs with ice cream than without it!" :
                                (Math.random() < 0.3?
                                "when you make ice cream with bugs you are saving the planet from starvation.":
                                        (Math.random() < 0.3?
                                 " love + best ingredients = BUGelato":
                                                "take a free BUGelato at any grocery store by showing your" +
                                                        " promotional code  I<3EATINGBUGS"+Math.round(1000*Math.random()))))+
                        "\n\n__________________________________________________________\n" +
                            "BUGelato Corp. Intellectual property - All Rights Reserved");
                break;
            }else{
                System.out.println("     It did not work well, please enter something I can understand.");
            }
        }

        }catch (IOException e){
            e.printStackTrace();
        }
    }


    private static double[][] fillA( double [][][][] x, int nprods, int ntimes, int nsales, int nproductions){
        double [][] ret = new double [ntimes][nsales];
            for(int j = 0; j<nsales; j++){
                for(int k = 0; k<ntimes; k++){
                        double val = 0;
                        for (int i = 0; i<nproductions; i++)
                            val += x[i][j][k][nprods-1];
                        ret[k][j] = val;
                }
            }
        return ret;
    }
    private static int[] computeAlpha(double [][][][] x,int  nsales, int nproductions,int  nprods, int ntimes) {

        int[] alpha = new int[nproductions];

        double[] suptot = new double[nsales];

        for (int i = 0; i < nproductions; i++) {
            for (int j = 0; j < nsales; j++) {
                for (int n = 0; n < nprods; n++) {
                    for (int k = 0; k < ntimes; k++) {
                        suptot[i] += x[i][j][k][n];
                    }
                }
            }

        }

        for (int i = 0; i < nproductions; i++) {
            if (suptot[i] != 0) {
                alpha[i] = 1;
            }
        }

        return alpha;
    }
    private static double computeDemand( double a, double b, double M, double P,double Q, double  A,double S, double p){
        double ret = 0;
        ret = M*(Math.pow(Math.max(2.0-p/P, 0.0), a))*(1.0+b*(0.5*S+A)/1e4)*Q/10;
        return ret;
    }
    private static int[]  computeDelta(double [][][][] x,int  nsales, int nproductions,int  nprods, int ntimes) {


        int[] delta = new int [nsales];
        double[] suptot = new double[nsales];

        for (int i = 0; i < nproductions; i++) {
            for (int j = 0; j < nsales; j++) {
                for (int n = 0; n < nprods; n++) {
                    for (int k = 0; k < ntimes; k++) {
                        suptot[j] += x[i][j][k][n];
                    }
                }
            }

        }

        for (int j = 0; j < nsales; j++) {
            if (suptot[j] != 0) {
                delta[j] = 1;
            }
        }
        return delta;

    }
    private static double[][] computeQuantities( double [][][][] x, double [][] conversion_factor, int  nsales, int nproductions,
                                                 int  nprods, int ntimes, int nrawmaterials){
        double [][] q = new double[nrawmaterials][nproductions];
        for(int m = 0; m<nrawmaterials; m++) {
            for (int i = 0; i < nproductions; i++) {
                for(int j = 0; j<nsales; j++){
                    for(int n = 0 ; n<nprods;n++){
                        for(int k = 0; k<ntimes; k++){
                            q[m][i] += x[i][j][k][n]*conversion_factor[m][n];
                        }
                    }
                }

            }
        }

        return q;
    }
    private static double[][] computeEnergies( double [][][][] x, double [][] energy_conversion, int  nsales, int nproductions,
                                                 int  nprods, int ntimes){
        double [][] energy_consumption = new double[nproductions][ntimes];

            for (int i = 0; i < nproductions; i++) {
                for(int j = 0; j<nsales; j++){
                    for(int n = 0 ; n<nprods;n++){
                        for(int k = 0; k<ntimes; k++){
                            energy_consumption[i][k] += x[i][j][k][n]*energy_conversion[n][i];
                        }
                    }
                }

            }


        return energy_consumption;
    }
    private static Object[] computeWharehouses(double [][][][] x, double [][][][] demand, double[][] p,
                                               int nmkts, int  nsales, int nproductions,int  nprods, int ntimes){
//todo verify consistency with respect to the center closed
        double [][][] arriving = new double [nprods][ntimes][nsales];
        double [][][] yinv = new double [nprods][ntimes][nsales];
        double [][][][] v = new double [nmkts][ntimes][nsales][nprods];



        for(int j = 0; j<nsales; j++){
            for(int k = 0; k<ntimes; k++){
                for(int n = 0 ; n<nprods; n++){
                    for(int i = 0; i<nproductions; i++){
                        arriving[n][k][j] += x[i][j][k][n];
                    }
                }
            }
        }

//k=0
        for(int j = 0; j<nsales; j++){

            for(int n = 0 ; n<nprods; n++){
                double availability = arriving[n][0][j];
                double demand0 = Math.floor(demand[0][j][0][n]);
                double demand1 = Math.floor(demand[1][j][0][n]);
                if(p[0][n] > p[1][n]){
                    double sold = Math.min( demand0, availability );
                    v[0][0][j][n] = sold;
                    double soldsecondround = Math.min( demand1, availability - sold);
                    v[1][0][j][n] = soldsecondround;
                    yinv[n][0][j] = availability - sold- soldsecondround;
                }else{
                    double sold = Math.min( demand1, availability );
                    v[1][0][j][n] = sold;
                    double soldsecondround = Math.min( demand0, availability - sold);
                    v[0][0][j][n] = soldsecondround;
                    yinv[n][0][j] = availability - sold- soldsecondround;
                }

            }
        }


        //k>0
        for (int k = 1; k<ntimes; k++) {
            for (int j = 0; j < nsales; j++) {

                for (int n = 0; n < nprods; n++) {
                    double availability = arriving[n][k][j];
                    double demand0 = Math.floor(demand[0][j][k][n]);
                    double demand1 = Math.floor(demand[1][j][k][n]);
                    double sold = 0;
                    double soldsecondround = 0;

                    if (p[0][n] > p[1][n]) {
                        sold = Math.min(demand0, availability);
                        v[0][k][j][n] = sold;
                        soldsecondround = Math.min(demand1, availability - sold);
                        v[1][k][j][n] = soldsecondround;
                    } else {
                        sold = Math.min(demand1, availability);
                        v[1][k][j][n] = sold;
                        soldsecondround = Math.min(demand0, availability - sold);
                        v[0][k][j][n] = soldsecondround;
                    }
                    yinv[n][k][j] = availability - sold - soldsecondround + yinv[n][k-1][j];
                }
            }
        }

        Object[] ret = new Object[2];
        ret[0] = yinv;
        ret[1] = v;
        return ret;
    }
    private static Object[] computeTrucks( double[][][] Ti, int ntrucks, int ntimes, int nproductions){
        double [][] Tsums = new double[ntimes][ntrucks];
        double [] Ttot = new double[ntrucks];
        double [] Treloc = new double[ntrucks];


        for(int t = 0; t<ntrucks; t++){
            for(int k = 0; k<ntimes; k++){
                for(int i = 0; i<nproductions; i++){
                    Tsums[k][t] += Ti[i][k][t];
                }
                if( Tsums[k][t] >= Ttot[t]){
                    Ttot[t] = Tsums[k][t];
                }
                if(k > 0 && (Tsums[k-1][t]-Tsums[k][t])>0){
                    Treloc[t] += (Tsums[k-1][t]-Tsums[k][t]);
                }
            }
        }



        Object[] ret = new Object[2];
        ret[0] = Ttot;
        ret[1] = Treloc;
        return ret;
    }
    private static double[] computeRevenues(double[][][][]v,double[][] p,int  nmkts,int  nsales,int  ntimes,int  nprods){

        double profit = 0;
        double profit0 = 0;
        double profit1 = 0;
        double [] ret = new double[3];

        for(int m = 0; m< nmkts; m++){
            for(int j = 0; j <nsales; j++){
                for(int k = 0; k<ntimes; k++){
                    for(int n = 0; n<nprods; n++){
                        profit += v[m][k][j][n]*p[m][n];
                        if (m==0){
                            profit0 += v[m][k][j][n]*p[m][n];
                        }else{
                            profit1 += v[m][k][j][n]*p[m][n];
                        }
                    }
                }
            }
        }

        ret[0] = profit;
        ret[1] = profit0;
        ret[2] = profit1;
        return ret;
    }
    private static double rawMaterialCosts(double[][]q, double[][]ar,double [][] br,int nvendors,int nrawmaterials,int nproductions){
        double cost = 0;
        double [][] auction = new double[nrawmaterials][nproductions];

        for(int r = 0; r<nrawmaterials; r++){
            for(int i = 0; i<nproductions; i++){
                auction[r][i] = Double.MAX_VALUE;
                for(int v= 0; v<nvendors; v++){
                    double offer = ar[r][v]*q[r][i]+br[r][v];
                    if( auction[r][i] > offer){
                        auction[r][i] = offer;
                    }
                }
                cost += auction[r][i];
            }
        }
        return cost;
    }
    private static double energyCosts(double[][] energy, double[]aenergy, double[] energy_threshold,
                                      int ntimes,int nproductions){
        double cost = 0;

        for(int k = 0; k<ntimes; k++){
            for(int i = 0; i<nproductions;i++){

                    double value;
                    if(energy[k][i] > energy_threshold[i]){
                        value = energy[i][k] - energy_threshold[i];
                        value = value*aenergy[1] + aenergy[0]*energy_threshold[i];
                    }else{
                        value = aenergy[0]*energy[i][k];
                    }

                cost += value;
            }
        }

        return cost;
    }
    private static double labourCost( double[][] wvar, double[] wfixed, double[] varpaychecks, double[] fixedpaychecks,
                                      int nproductions, int ntimes){
        double cost = 0;
        for(int i = 0; i<nproductions; i++){
            cost += wfixed[i]*fixedpaychecks[i];
            for(int k = 0; k<ntimes; k++){
                cost += wvar[k][i]*varpaychecks[i];
            }
        }

        return cost;
    }
    private static double fixedSales(int [] delta, double [] fixedsales, int nsales){
        double cost =0;
        for (int j = 0; j<nsales; j++){
            cost += delta[j]*fixedsales[j];
        }
        return cost;
    }
    private static double fixedPruductions(int [] alpha, double [] fixedproductions, int nproductions){
        double cost =0;
        for (int j = 0; j<nproductions; j++){
            cost += alpha[j]*fixedproductions[j];
        }
        return cost;
    }
    private static double transportationCost( double[][][][] x, double[][] transp_cost, double [] volume,
                                              int nproductions,int nsales,int ntimes,int nprods){
        double cost = 0;
        for(int i = 0; i<nproductions; i++)
            for(int j = 0; j<nsales; j++)
                for(int k = 0; k<ntimes; k++)
                    for(int n= 0; n<nprods; n++)
                        cost += x[i][j][k][n]*volume[n]*transp_cost[i][j];
        return cost;
    }
    private static double inventoryCost(double [][][] yinv, double[][] inventory_cost_param, int nsales, int ntimes, int nprods){
        double cost = 0;

        for(int j = 0; j<nsales; j++){
            for(int k = 0; k<ntimes; k++)
                for(int n = 0; n<nprods; n++)
                    cost += yinv[n][k][j]*inventory_cost_param[j][n];
        }

        return cost;
    }
    private static double fixedTruckCost( double [] Ttot, double[] truck_cost, int ntrucks){
        double cost = 0;

        for(int t = 0; t<ntrucks; t++){
            cost += Ttot[t]*truck_cost[t];
        }
        return cost;
    }
    private static double relocationCost( double [] Treloc, double[]  reloc_cost,int  ntrucks){
        double cost = 0;

        for(int t = 0; t<ntrucks; t++)
            cost += Treloc[t]*reloc_cost[t];

        return cost;

    }
    private static  double activeTruckCost( double [][][][] T, double [][] activation_truck_cost,
                                            int nproductions, int nsales, int ntimes, int ntrucks ){
        double cost = 0;
        for(int i = 0; i<nproductions; i++)
            for(int j = 0; j<nsales; j++)
                for(int k = 0; k<ntimes; k++)
                    for(int t= 0 ; t<ntrucks; t++)
                        cost += T[i][j][k][t]*activation_truck_cost[i][t];

        return cost;

    }
    private static double Sum(double[][] inv){
        double cost = 0;
        for(int i = 0 ; i<inv.length; i++)
            for(int j= 0; j<inv[0].length; j++)
                cost +=inv[i][j];

        return cost;
    }

    private static Object[] checkSales(int [] delta){
        int sum = 0;
        for( int i= 0; i<delta.length ; i++)
            sum += delta[i];
        Object[] ret = new Object[2];
        if (sum <= 7){
            ret[0] = true;
            return ret;
        }else{
            ret[0] = false;
            String msg = "INFEASIBILITY CHECK: all the Sales centers are open!";
            ret[1] = msg;
            return ret;
        }
    }
    private static Object[] checkNonNegativity( double [][][][]x,double[][][][] T,double[][][] Ti,
                                                double[] wfixed, double[][] wvar, double [][] p,
                                                double [][]investment,
                                                int nproductions, int nsales, int ntimes, int nprods, int nmkts,int ntrucks){




        Object[] ret = new Object[2];
        boolean check = true;
        String s = "";

        megaloop:
        for(int i = 0; i<nproductions; i++){
            if(wfixed[i] <= -1e-6){
                check = false;
                s = " fixed workers ";
                break megaloop;
            }
            for(int k = 0; k<ntimes; k++){
                if(wvar[k][i] <= -1e-6){
                    check = false;
                    s = " variable workers ";
                    break megaloop;
                }

                for(int t = 0; t<ntrucks; t++){
                    if(Ti[i][k][t] <= -1e-6){
                        check = false;
                        s = " allocated trucks ";
                        break megaloop;
                    }
                }
            }


            for(int j = 0; j<nsales; j++){
                for(int k = 0; k<ntimes; k++){
                    for(int n=0; n<nprods; n++){
                        if(x[i][j][k][n] <= -1e-6){
                            check = false;
                            s = " products ";
                            break megaloop;
                        }
                    }

                    for(int t = 0 ; t<ntrucks; t++){
                        if(T[i][j][k][t] <= -1e-6){
                            check = false;
                            s = " scheduled trucks ";
                            break megaloop;
                        }
                    }



                }
            }
        }


        miniloop:
        for(int m = 0; m<nmkts; m++){
            for(int n = 0; n<nprods; n++){
                if(p[m][n] <= -1e-6){
                    check = false;

                    s = " prices ";
                    break miniloop;
                }
            }
            for(int k = 0; k<ntimes; k++){
                if(investment[m][k] <= -1e-6){
                    check = false;
                    s = " investments ";
                    break miniloop;
                }
            }

        }

        ret[0] = check;
        if (check){
            return ret;
        }else{
            String msg = "INFEASIBILITY CHECK: some of the input variables is negative. Please check the following input: "+s+".";
            ret[1] = msg;
            return ret;
        }

    }
    private static  Object[] checkWorkersAndSaleCenters( int [] alpha, double[][][][] x,double [] wfixed,double[][] wvar,
                                                         double [][] worktime, double[] lb,double[] ub,
                                                         double hpermonth, int nproductions,int  nsales,
                                                         int ntimes, int nprods){
        double [][] supalpha = new double[ntimes][nproductions];
        Object[] ret = new Object[2];
        boolean checkbool = true;
        String s = "";

        loop:
        for(int k = 0; k<ntimes; k++){

            for(int i = 0; i<nproductions; i++){
                double lhs = 0;
                for(int j = 0; j<nsales; j++){
                    for(int n=0; n<nprods; n++){
                        lhs += x[i][j][k][n]*worktime[n][i];
                    }
                }
                supalpha[k][i] = lhs;


                double rhs = hpermonth*(wfixed[i]+wvar[k][i]);
                if ((rhs +1e-6 )<= lhs ){
                    checkbool = false;
                    s = " too few workers to cover the production requirements.";
                    break loop;
                }

            }
        }

        if(!checkbool){
            ret[0] = checkbool;
            s = "INFEASIBILITY CHECK:"+s;
            ret[1] = s;
            return ret;
        }


        loop:
        for(int i =0; i<nproductions; i++){
            for(int k = 0; k<ntimes; k++){
                double rhs = (wfixed[i]+wvar[k][i]);
                if( rhs >= (ub[i]*alpha[i] + 1e-6) || rhs <= lb[i]*alpha[i]-1e-6){
                    checkbool = false;
                    s = " the solution exceeds the bounds.";
                    break loop;
                }
            }
        }


        if(!checkbool){
            ret[0] = checkbool;
            s = "INFEASIBILITY CHECK:"+s;
            ret[1] = s;
            return ret;
        }else{
            ret[0] = checkbool;
            return ret;
        }
    }
    private static Object[] checkTrucks( double [][][][]x, double [][][][]T,double [] truck_volume,double [] volume,
                                         int nproductions, int nsales,  int ntimes,  int nprods,int ntrucks){
        Object[] ret = new Object[2];
        String s="";
        boolean checkbool = true;
        ret[0] = checkbool;




        for(int k = 0; k<ntimes; k++) {
            for (int i = 0; i < nproductions; i++) {
                for (int j = 0; j < nsales; j++) {
                    double lhs = 0;
                    for (int n = 0; n < nprods; n++) {
                        lhs += volume[n]*x[i][j][k][n];
                    }
                    double rhs = 0;
                    for (int t = 0; t<ntrucks; t++){
                        rhs += T[i][j][k][t]*truck_volume[t];
                    }

                    if( lhs >= (rhs + 1e-6)){
                        checkbool = false;
                        s = "INFEASIBILITY CHECK: the trucks provided do not cover the production load.";
                        ret[0] = checkbool;
                        ret[1] = s;
                        return ret;
                    }
                }
            }
        }

        return ret;
    }

    private static Object[] fillInput(int nproductions, int nsales, int ntimes, int nprods, int nmkts, int ntrucks, String inputpath){
        double [][][][] x = new double [nproductions][nsales][ntimes][nprods];
        double [][] p = new double [nmkts][nprods];
        double [][] wvar = new double [ntimes][nproductions];
        double [] wfixed = new double [nproductions];
        double [][] investment = new double [nmkts][ntimes];
        double [][][][] T = new double[nproductions][nsales][ntimes][ntrucks];
        double [][][] Ti = new double[nproductions][ntimes][ntrucks];



        Object[] ret = new Object[7];
        ret[0] = x;
        ret[1] = p;
        ret[2] = wvar;
        ret[3] = wfixed;
        ret[4] = investment;
        ret[5] = T;
        ret[6] = Ti;
        return ret;
    }
    private static Object[] fillParameters(int nproductions, int nsales, int ntimes, int nprods, int nmkts,
                                           int ntrucks, int nrawmaterials, int nvendors){
        double hpermonth = 20*8;
        double [][] worktime = new double [nprods][nproductions];
        double [] fixedpaycheck = new double[]{1500, 1200, 1250};
        double [] variablepaycheck = new double[]{1900, 1800, 1850};
        double [] lb ={15, 20, 15};// new double [nproductions];
        double [] ub ={30, 30, 35};

        double [] volume = new double [] {.6, .3, .8, .4, .01};
        double [][] transp_cost = new double [nproductions][nsales];
        double [] reloc_cost = new double [ntrucks];
        double [] truck_cost = new double [ntrucks];
        double [][] activation_truck_cost = new double [nproductions][ntrucks];
        double [] truck_volume = new double[ntrucks];

        double [][] conversion_factor = new double [nrawmaterials][nprods];
        double [][] ar = new double[nrawmaterials][nvendors];
        double [][] br = new double[nrawmaterials][nvendors];

        double [][] energy_conversion = new double[nprods][nproductions];
        double [] aenergy = new double[] {.1, .2};
        double [] energy_threshold = new double[] {56286.7, 6940.54, 8100.66};

        double [] fixedsales = new double[]{20000, 18000, 10000, 9000, 12000, 15000, 22000, 13000};
        double [] fixedproductions = new double[] {800000, 600000, 650000};

        double[][] inventory_cost_param = new double[nsales][nprods];


        //demand params
        double [][] M = new double[nsales][nprods];
        double [][] P = new double[nsales][nprods];
        double a = 0;
        double [][] Q = new double[nmkts][ntimes];
        double b = 0;




        Object [] ret = new Object[27];
        ret[0] = hpermonth;
        ret[1] = worktime;
        ret[2] = fixedpaycheck;
        ret[3] = variablepaycheck;
        ret[4] = lb;
        ret[5] = ub;
        ret[6] = volume;
        ret[7] = transp_cost;
        ret[8] = reloc_cost;
        ret[9] = truck_cost;
        ret[10] = activation_truck_cost;
        ret[11] = truck_volume;
        ret[12] = conversion_factor;
        ret[13] = ar;
        ret[14] = br;
        ret[15] = energy_conversion;
        ret[16] = aenergy;
        ret[17] = energy_threshold;
        ret[18] = fixedsales;
        ret[19] = fixedproductions;
        ret[20] = inventory_cost_param;
        ret[21] = M;
        ret[22] = P;
        ret[23] = a;
        ret[24] = Q;
        ret[25] = b;
        return ret;
    }
}
