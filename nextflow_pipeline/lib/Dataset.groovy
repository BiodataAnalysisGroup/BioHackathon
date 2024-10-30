// Dataset class definition
class Dataset {
    String name
	String type
	String path

    Dataset(String name, String type, String path) {
        this.name 	= name
        this.type  	= type
		this.path	= path
    }

    String toString() {
        "$name"
    }
}