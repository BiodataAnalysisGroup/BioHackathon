// Perturbation task class definition
class Task {
    String name
    List<String> compatibleMetrics

    Task(String name, List<String> compatibleMetrics) {
        this.name = name
        this.compatibleMetrics = compatibleMetrics
    }

    // Check if the metric is compatible with a given perturbation task
    boolean isMetricCompatible(String metric) {
        return this.compatibleMetrics.contains(metric)
    }

    String toString() {
        "$name"
    }
}