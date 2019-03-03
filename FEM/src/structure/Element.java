package structure;

public class Element {
    private int id;
    public int[] nodes;

    public double[][] hMatrix;
    public double[][] cMatrix;
    public double[] pVector;

    private static final double[] KSI = {-1 / Math.sqrt(3), 1 / Math.sqrt(3), 1 / Math.sqrt(3), -1 / Math.sqrt(3)};
    private static final double[] ETA = {-1 / Math.sqrt(3), -1 / Math.sqrt(3), 1 / Math.sqrt(3), 1 / Math.sqrt(3)};
    private static double[][] SHAPE_FUNCTION;
    private static double[][] dN_dKSI;
    private static double[][] dN_dETA;

    public Element() {
        nodes = new int[4];
        hMatrix = new double[4][4];
        cMatrix = new double[4][4];
        pVector = new double[4];


        if (SHAPE_FUNCTION == null) {

            SHAPE_FUNCTION = new double[4][4];
            dN_dKSI = new double[4][4];
            dN_dETA = new double[4][4];

            for (int i = 0; i < 4; ++i) {
                SHAPE_FUNCTION[0][i] = 0.25 * (1 - KSI[i]) * (1 - ETA[i]);
                SHAPE_FUNCTION[1][i] = 0.25 * (1 + KSI[i]) * (1 - ETA[i]);
                SHAPE_FUNCTION[2][i] = 0.25 * (1 + KSI[i]) * (1 + ETA[i]);
                SHAPE_FUNCTION[3][i] = 0.25 * (1 - KSI[i]) * (1 + ETA[i]);
            }
            for (int i = 0; i < 4; ++i) {
                dN_dKSI[0][i] = -0.25 * (1 - ETA[i]);
                dN_dKSI[1][i] = 0.25 * (1 - ETA[i]);
                dN_dKSI[2][i] = 0.25 * (1 + ETA[i]);
                dN_dKSI[3][i] = -0.25 * (1 + ETA[i]);
            }
            for (int i = 0; i < 4; ++i) {
                dN_dETA[0][i] = -0.25 * (1 - KSI[i]);
                dN_dETA[1][i] = -0.25 * (1 + KSI[i]);
                dN_dETA[2][i] = 0.25 * (1 + KSI[i]);
                dN_dETA[3][i] = 0.25 * (1 - KSI[i]);
            }
        }
    }
    public void countElement(Node one, Node two, Node three, Node four, boolean[] bc, double k, double c, double ro, double alfa, double tEnv) {
        double[][] jMatrix = new double[4][4];
        double[] detJ = new double[4];
        double[][] dn_dx = new double[4][4];
        double[][] dn_dy = new double[4][4];

        for (int pc = 0; pc < 4; ++pc) {
            jMatrix[0][pc] = dN_dETA[0][pc] * one.getY() + dN_dETA[1][pc] * two.getY() + dN_dETA[2][pc] * three.getY() + dN_dETA[3][pc] * four.getY();//dy/dn
            jMatrix[1][pc] = -(dN_dKSI[0][pc] * one.getY() + dN_dKSI[1][pc] * two.getY() + dN_dKSI[2][pc] * three.getY() + dN_dKSI[3][pc] * four.getY());//-dy/de
            jMatrix[2][pc] = -(dN_dETA[0][pc] * one.getX() + dN_dETA[1][pc] * two.getX() + dN_dETA[2][pc] * three.getX() + dN_dETA[3][pc] * four.getX());//-dx/dn
            jMatrix[3][pc] = dN_dKSI[0][pc] * one.getX() + dN_dKSI[1][pc] * two.getX() + dN_dKSI[2][pc] * three.getX() + dN_dKSI[3][pc] * four.getX();//dx/de
        }
        for (int pc = 0; pc < 4; ++pc) {
            detJ[pc] = jMatrix[0][pc] * jMatrix[3][pc] - jMatrix[1][pc] * jMatrix[2][pc];
        }
        for (int pc = 0; pc < 4; ++pc) {
            for (int i = 0; i < 4; ++i) {
                dn_dx[i][pc] = jMatrix[0][pc] / detJ[pc] * dN_dKSI[i][pc] + jMatrix[1][pc] / detJ[pc] * dN_dETA[i][pc];
                dn_dy[i][pc] = jMatrix[2][pc] / detJ[pc] * dN_dKSI[i][pc] + jMatrix[3][pc] / detJ[pc] * dN_dETA[i][pc];
            }
        }

        for (int pc = 0; pc < 4; ++pc) {
            for (int i = 0; i < 4; ++i) {
                for (int j = 0; j < 4; ++j) {
                    hMatrix[i][j] += k * (dn_dx[i][pc] * dn_dx[j][pc] * detJ[pc] + dn_dy[i][pc] * dn_dy[j][pc] * detJ[pc]);
                }
            }
        }

        for (int pc = 0; pc < 4; ++pc) {
            for (int i = 0; i < 4; ++i) {
                for (int j = 0; j < 4; ++j) {
                    cMatrix[i][j] += SHAPE_FUNCTION[i][pc] * SHAPE_FUNCTION[j][pc] * detJ[pc] * c * ro;
                }
            }
        }

        if (bc[0] == true) {
            double L = two.getX() - one.getX();
            hMatrix[0][0] += alfa * L / 2 * (0.25 * (1 + 1 / Math.sqrt(3)) * (2));//ONLY FOR RECTANGULAR ELEMENT
            hMatrix[0][1] += alfa * L / 2 * (0.25 * (1 - 1 / Math.sqrt(3)) * (2));
            hMatrix[1][0] += alfa * L / 2 * (0.25 * (1 + 1 / Math.sqrt(3)) * (2));
            hMatrix[1][1] += alfa * L / 2 * (0.25 * (1 - 1 / Math.sqrt(3)) * (2));
            pVector[0] += alfa * tEnv * L / 2 * (shapeFunction1(-1 / Math.sqrt(3), -1) + shapeFunction1(1 / Math.sqrt(3), -1));
            pVector[1] += alfa * tEnv * L / 2 * (shapeFunction2(-1 / Math.sqrt(3), -1) + shapeFunction2(1 / Math.sqrt(3), -1));
        }

        if (bc[1] == true) {
            double L = three.getY() - two.getY();
            hMatrix[1][1] += alfa * L / 2 * (0.25 * (2) * (1 + 1 / Math.sqrt(3)));//ONLY FOR RECTANGULAR ELEMENT
            hMatrix[1][2] += alfa * L / 2 * (0.25 * (2) * (1 - 1 / Math.sqrt(3)));
            hMatrix[2][1] += alfa * L / 2 * (0.25 * (2) * (1 - 1 / Math.sqrt(3)));
            hMatrix[2][2] += alfa * L / 2 * (0.25 * (2) * (1 + 1 / Math.sqrt(3)));
            pVector[1] += alfa * tEnv * L / 2 * (shapeFunction2(1, -1 / Math.sqrt(3)) + shapeFunction2(1, 1 / Math.sqrt(3)));
            pVector[2] += alfa * tEnv * L / 2 * (shapeFunction3(1, -1 / Math.sqrt(3)) + shapeFunction3(1, 1 / Math.sqrt(3)));
        }

        if (bc[2] == true) {
            double L = three.getX() - four.getX();
            hMatrix[2][2] += alfa * L / 2 * (0.25 * (1 + 1 / Math.sqrt(3)) * (2));//ONLY FOR RECTANGULAR ELEMENT
            hMatrix[2][3] += alfa * L / 2 * (0.25 * (1 - 1 / Math.sqrt(3)) * (2));
            hMatrix[3][2] += alfa * L / 2 * (0.25 * (1 - 1 / Math.sqrt(3)) * (2));
            hMatrix[3][3] += alfa * L / 2 * (0.25 * (1 + 1 / Math.sqrt(3)) * (2));
            pVector[2] += alfa * tEnv * L / 2 * (shapeFunction3(1 / Math.sqrt(3), 1) + shapeFunction3(-1 / Math.sqrt(3), 1));
            pVector[3] += alfa * tEnv * L / 2 * (shapeFunction4(1 / Math.sqrt(3), 1) + shapeFunction4(-1 / Math.sqrt(3), 1));
        }

        if (bc[3] == true) {
            double L = three.getY() - two.getY();
            hMatrix[0][0] += alfa * L / 2 * (0.25 * (2) * (1 + 1 / Math.sqrt(3)));//ONLY FOR RECTANGULAR ELEMENT
            hMatrix[0][3] += alfa * L / 2 * (0.25 * (2) * (1 - 1 / Math.sqrt(3)));
            hMatrix[3][0] += alfa * L / 2 * (0.25 * (2) * (1 - 1 / Math.sqrt(3)));
            hMatrix[3][3] += alfa * L / 2 * (0.25 * (2) * (1 + 1 / Math.sqrt(3)));
            pVector[0] += alfa * tEnv * L / 2 * (shapeFunction1(-1, 1 / Math.sqrt(3)) + shapeFunction1(-1, -1 / Math.sqrt(3)));
            pVector[3] += alfa * tEnv * L / 2 * (shapeFunction4(-1, 1 / Math.sqrt(3)) + shapeFunction4(-1, -1 / Math.sqrt(3)));
        }
    }
    public void showLocalH() {
        System.out.println("Local H matrix for " + this.id + "element:");
        for (int i = 0; i < 4; ++i) {
            for (int j = 0; j < 4; ++j) {
                System.out.print(hMatrix[i][j] + " ");
            }
            System.out.println();
        }
        System.out.println();
    }
    public void showLocalC() {
        System.out.println("Local C matrix for " + this.id + "element:");
        System.out.println();
        for (int i = 0; i < 4; ++i) {
            for (int j = 0; j < 4; ++j) {
                System.out.print(cMatrix[i][j] + " ");
            }
            System.out.println();
        }
        System.out.println();
    }
    public void showLocalP() {
        System.out.println("Local P vector for " + this.id + "element:");
        System.out.println();
        for (int i = 0; i < 4; ++i) {
            System.out.println(pVector[i]);
        }
        System.out.println();
    }
    public void setId(int id) {
        this.id = id;
    }
    public void setNodes(int one, int two, int three, int four) {
        this.nodes[0] = one;
        this.nodes[1] = two;
        this.nodes[2] = three;
        this.nodes[3] = four;
    }
    private static double shapeFunction1(double ksi, double eta) {
        return (0.25 * (1 - ksi) * (1 - eta));
    }
    private static double shapeFunction2(double ksi, double eta) {
        return (0.25 * (1 + ksi) * (1 - eta));
    }
    private static double shapeFunction3(double ksi, double eta) {
        return (0.25 * (1 + ksi) * (1 + eta));
    }
    private static double shapeFunction4(double ksi, double eta) {
        return (0.25 * (1 - ksi) * (1 + eta));
    }

}
