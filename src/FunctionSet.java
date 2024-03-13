import java.util.function.Function;

class FunctionSet {
    public static double weierstrass_function(double x) {
        int n = 5;
        double b = 0.5;
        int a = 13;

        double f_x = 0;
        for (int i = 0; i < n; i++) {
            f_x += Math.pow(b, i) * Math.cos(Math.pow(a, n) * Math.PI * x);
        }

        return f_x;
    }

    public static double gamma_function(double x) {
        double tmp = (x - 0.5) * Math.log(x + 4.5) - (x + 4.5);
        double ser = 1.0 +
                76.18009173 / (x + 0.0) - 86.50532033 / (x + 1.0) +
                24.01409822 / (x + 2.0) - 1.231739516 / (x + 3.0) +
                0.00120858003 / (x + 4.0) - 0.00000536382 / (x + 5.0);
        return Math.exp(tmp + Math.log(ser * Math.sqrt(2 * Math.PI)));
    }

    public static Function<Double, Double> get_function(int n) {
        return switch (n) {
            case (1) -> FunctionSet::weierstrass_function;
            case (2) -> FunctionSet::gamma_function;
            default -> throw new UnsupportedOperationException("Function " + n + " not defined.");
        };
    }
}
