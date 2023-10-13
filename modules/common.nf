import fr.ens.biologie.genomique.kenetre.io.CompressionType
import fr.ens.biologie.genomique.kenetre.storage.FileStorage
import fr.ens.biologie.genomique.kenetre.storage.FileGenomeDescStorage
import fr.ens.biologie.genomique.kenetre.bio.GenomeDescription
import fr.ens.biologie.genomique.kenetre.log.DummyLogger
import fr.ens.biologie.genomique.kenetre.log.StandardErrorLogger

// The logger to use
//logger = new StandardErrorLogger("stderr", [:])
logger = new DummyLogger()

def read_conf() {

    result = [:]
    conf_path = System.getProperty("user.home") + "/.eoulsan-nf"

    File conf_file = new File(conf_path)

    if (!conf_file.exists()) {
        return result
    }

    content = conf_file.eachLine {
        line -> {
            fields = line.split('=')
            if (fields.length ==2 ) {
            result[fields[0].trim()] = fields[1].trim()
            }
        }
    }
    return result
}

def get_path(p, storages) {

    if (p.toLowerCase().startsWith("genome://") && storages.containsKey("main.genome.storage.path") ) {
        fs = new FileStorage(storages["main.genome.storage.path"], [".fasta", ".fa", ".fna"])
        return fs.getFile(p)
    } else if (p.toLowerCase().startsWith("gff://") && storages.containsKey("main.gff.storage.path") ) {
        fs = new FileStorage(storages["main.gff.storage.path"], [".gff", ".gff3"])
        return fs.getFile(p)
    } else if (p.toLowerCase().startsWith("gtf://") && storages.containsKey("main.gtf.storage.path") ) {
        fs = new FileStorage(storages["main.gtf.storage.path"], [".gtf"])
        return fs.getFile(p)
    } else if (p.toLowerCase().startsWith("additionalannotation://") && storages.containsKey("main.additional.annotation.storage.path") ) {
        fs = new FileStorage(storages["main.additional.annotation.storage.path"], [".tsv", "txt"])
        return fs.getFile(p)
    }

    return new File(p)
}

def get_genome_desc(p, storages) {

    if (storages == null) {
        return p
    }

    unix_file = get_path(p, storages)

    if (storages.containsKey("main.genome.desc.storage.path") ) {
        gds = FileGenomeDescStorage.getInstance(storages["main.genome.desc.storage.path"], logger)

        result = gds.get(unix_file.toString())
        if (result != null) {
            return result
        }
    }

    // TODO Store new genome description in repository
    return GenomeDescription.createGenomeDescFromFasta(CompressionType.open(unix_file), unix_file.toString())
}

def input_stream(p) {

    ct = CompressionType.getCompressionTypeByFile(p)
    return ct.open(p)
}

def output_stream(p) {

    ct = CompressionType.getCompressionTypeByFile(p)
    return ct.create(p)
}