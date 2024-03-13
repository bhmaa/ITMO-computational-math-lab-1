import java.io.*;

class Main {
    public static void main(String[] args) throws IOException {
        BufferedReader bufferedReader = new BufferedReader(new InputStreamReader(System.in));
        BufferedWriter bufferedWriter = new BufferedWriter(new OutputStreamWriter(System.out));

        int f = Integer.parseInt(bufferedReader.readLine().trim());
        double a = Double.parseDouble(bufferedReader.readLine().trim());
        double b = Double.parseDouble(bufferedReader.readLine().trim());
        double x = Double.parseDouble(bufferedReader.readLine().trim());

        double result = Result.interpolate_by_spline(f, a, b, x);

        bufferedWriter.write(String.valueOf(result));
        bufferedWriter.newLine();

        bufferedReader.close();
        bufferedWriter.close();
    }
}
