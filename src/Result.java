import java.util.ArrayList;
import java.util.Collections;
import java.util.List;
import java.util.function.Function;

class Result {
    private static List<Double> get_chebyshev_nodes(double start_point, double end_point, int nodes_count) {
        List<Double> nodes = new ArrayList<>();
        for (int k = 1; k <= nodes_count; k++) {
            nodes.add(0.5 * (start_point + end_point) + 0.5 * (end_point - start_point) * Math.cos((2.0 * k - 1.0)
                    / (2.0 * nodes_count) * Math.PI));
        }
        Collections.reverse(nodes);
        return nodes;
    }

    private static class CubicSplineCoefficients {
        double a, b, c, d;

        public CubicSplineCoefficients(double a, double b, double c, double d) {
            this.a = a;
            this.b = b;
            this.c = c;
            this.d = d;
        }
    }

    private static CubicSplineCoefficients[] calculate_cubic_splines_coefficients(List<Double> x, List<Double> y) {
        int intervals_count = x.size() - 1;
        double[] lower_matrix_diagonal = new double[intervals_count];
        lower_matrix_diagonal[0] = 1;
        double[] upper_matrix_diagonal = new double[intervals_count - 1];
        upper_matrix_diagonal[0] = 0;

        for (int i = 1; i < intervals_count; ++i) {
            lower_matrix_diagonal[i] = (2 * (x.get(i + 1) - x.get(i - 1))) - (x.get(i) - x.get(i - 1))
                    * upper_matrix_diagonal[i - 1];
            if (i != intervals_count - 1) {
                upper_matrix_diagonal[i] = (x.get(i + 1) - x.get(i)) / lower_matrix_diagonal[i];
            }
        }
        double[] forward_solution_for_cs = new double[intervals_count];
        forward_solution_for_cs[0] = 0;
        for (int i = 1; i < intervals_count; ++i) {
            forward_solution_for_cs[i] =
                    ((3 * (y.get(i + 1) - y.get(i)) / (x.get(i + 1) - x.get(i)) - 3 * (y.get(i) - y.get(i - 1))
                            / (x.get(i) - x.get(i - 1))) - (x.get(i) - x.get(i - 1)) * upper_matrix_diagonal[i - 1])
                            / lower_matrix_diagonal[i];
        }
        CubicSplineCoefficients[] coefficients = new CubicSplineCoefficients[intervals_count];
        for (int i = 0; i < intervals_count; i++) {
            coefficients[i] = new CubicSplineCoefficients(y.get(i), 0, forward_solution_for_cs[intervals_count - 1], 0);
        }

        for (int j = intervals_count - 2; j >= 0; --j) {
            coefficients[j].c = forward_solution_for_cs[j] - upper_matrix_diagonal[j] * coefficients[j + 1].c;
            coefficients[j].b = (y.get(j + 1) - y.get(j)) / (x.get(j + 1) - x.get(j)) - (x.get(j + 1) - x.get(j))
                    * (2 * coefficients[j].c + coefficients[j + 1].c) / 3;
            coefficients[j].d = (coefficients[j + 1].c - coefficients[j].c) / (3 * (x.get(j + 1) - x.get(j)));
        }
        return coefficients;
    }

    private static double get_cubic_spline_value(List<CubicSplineCoefficients> coefficients, List<Double> x, double x_value) {
        if (x_value > x.get(coefficients.size()) || x_value < x.get(0)) {
            throw new RuntimeException("Точка не является промежуточным значением для полученного датасета. "
                    + "Пределы полученного датасета: [" + x.get(0) + ", " + x.get(coefficients.size()) + "]");
        }
        int interval = 0;
        while (interval < coefficients.size() && x_value > x.get(interval + 1)) {
            ++interval;
        }
        double delta_x = x_value - x.get(interval);
        return coefficients.get(interval).a +
                coefficients.get(interval).b * delta_x +
                coefficients.get(interval).c * Math.pow(delta_x, 2) +
                coefficients.get(interval).d * Math.pow(delta_x, 3);
    }

    public static double interpolate_by_spline(int f, double a, double b, double x) {
        Function<Double, Double> interpolated_function = FunctionSet.get_function(f);
        int nodes_count = 5;
        List<Double> chebyshev_nodes = get_chebyshev_nodes(a, b, nodes_count);
        List<Double> function_values = new ArrayList<>();
        for (double node : chebyshev_nodes) {
            function_values.add(interpolated_function.apply(node));
        }
        List<CubicSplineCoefficients> splines = List.of(calculate_cubic_splines_coefficients(chebyshev_nodes,
                function_values));
        double interpolated_value = get_cubic_spline_value(splines, chebyshev_nodes, x);
        if (Double.isNaN(interpolated_value)) {
            throw new RuntimeException("Похоже, что это нельзя решить...");
        }
        double function_value = interpolated_function.apply(x);
        while (Math.abs(function_value - interpolated_value) >= 0.01) {
            chebyshev_nodes.clear();
            function_values.clear();
            nodes_count = nodes_count + 5;
            chebyshev_nodes = get_chebyshev_nodes(a, b, nodes_count);
            for (double node : chebyshev_nodes) {
                function_values.add(interpolated_function.apply(node));
            }
            splines = List.of(calculate_cubic_splines_coefficients(chebyshev_nodes, function_values));
            interpolated_value = get_cubic_spline_value(splines, chebyshev_nodes, x);
            if (Double.isNaN(interpolated_value) || nodes_count > 10000) {
                throw new RuntimeException("Похоже, что это нельзя решить...");
            }
        }
        return interpolated_value;
    }
}
