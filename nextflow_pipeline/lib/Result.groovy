// Result class definition
class Result {
    Run run    
    String metric

    Result(Run run, String metric){
        this.run    = run
        this.metric = metric
    }

    String toString() {
        "Evaluation of tool $run.name on dataset $run.dataset using metric $metric"
    }
}