import java.util.*;
import java.io.*;

public class LinearEquationSolver {
   public static void main(String[] args) {
      Scanner scanner = new Scanner(System.in);
      
      System.out.print("Enter the number of linear equations (n <= 10): ");
      int n = scanner.nextInt();
      
      double[][] A = new double[n][n];
      double[] b = new double[n];
      
      System.out.print("Enter coefficients by 'command' or from a 'file': ");
      String inputType = scanner.next();
      
      if (inputType.equalsIgnoreCase("command")) {
         System.out.println("Enter each row with coefficients followed by the constant term:");
         for (int i = 0; i < n; i++) {
               for (int j = 0; j < n; j++) {
                  A[i][j] = scanner.nextDouble();
               }
               b[i] = scanner.nextDouble();
         }
      } else {
         System.out.print("Enter the file name with the coefficients: ");
         String fileName = scanner.next();
         Scanner fileScanner = null;
         try {
               fileScanner = new Scanner(new File(fileName));
               for (int i = 0; i < n; i++) {
                  for (int j = 0; j < n; j++) {
                     A[i][j] = fileScanner.nextDouble();
                  }
                  b[i] = fileScanner.nextDouble();
               }
         } catch (FileNotFoundException e) {
               System.out.println("File not found.");
               return;
         } finally {
               if (fileScanner != null) {
                  fileScanner.close();
               }
         }
      }

      System.out.print("Enter the desired stopping error: ");
      double tolerance = scanner.nextDouble();
      
      double[] xInit = new double[n];
      System.out.print("Enter the starting solution (e.g., '0 0 0'): ");
      for (int i = 0; i < n; i++) {
         xInit[i] = scanner.nextDouble();
      }

      System.out.println("\nJacobi Method:");
      double[] jacobiResult = iterativeMethod(A, b, xInit, tolerance, true);
      System.out.println("Final result from Jacobi Method: " + Arrays.toString(jacobiResult));
      
      System.out.println("\nGauss-Seidel Method:");
      double[] gaussSeidelResult = iterativeMethod(A, b, xInit, tolerance, false);
      System.out.println("Final result from Gauss-Seidel Method: " + Arrays.toString(gaussSeidelResult));
      
      scanner.close();
   }
   
   private static double[] iterativeMethod(double[][] A, double[] b, double[] xInit, double tolerance, boolean isJacobi) {
      int n = A.length;
      double[] x = Arrays.copyOf(xInit, n);
      double[] xNew = new double[n];
      boolean converged = false;
      int iteration = 0;
      
      while (!converged && iteration < 50) {
         for (int i = 0; i < n; i++) {
               double sum = 0.0;
               for (int j = 0; j < n; j++) {
                  if (j != i) {
                     sum += A[i][j] * (isJacobi ? xInit[j] : x[j]);
                  }
               }
               xNew[i] = (b[i] - sum) / A[i][i];
         }
         
         double error = computeError(x, xNew);
         System.out.println("Iteration " + (iteration + 1) + ": " + Arrays.toString(xNew) + "    " + error);
         
         if (error < tolerance) {
               converged = true;
         }
         
         System.arraycopy(xNew, 0, x, 0, n);
         if (isJacobi) {
               System.arraycopy(xNew, 0, xInit, 0, n);
         }
         
         iteration++;
      }
      
      if (!converged) {
         System.out.println("Maximum iterations reached without achieving the desired accuracy.");
      }
      return x;
   }
   
   private static double computeError(double[] xOld, double[] xNew) {
      double norm = 0.0;
      for (int i = 0; i < xOld.length; i++) {
         norm += (xNew[i] - xOld[i]) * (xNew[i] - xOld[i]);
      }
      return Math.sqrt(norm);
   }
}