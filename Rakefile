# Copyright 2015 Anurag Priyam - MIT License
#
# Same species, annotation lift over pipeline.
#
# Based on the lift over procedure deseribed at:
# http://genomewiki.ucsc.edu/index.php/LiftOver_Howto &
# http://hgwdev.cse.ucsc.edu/~kent/src/unzipped/hg/doc/liftOver.txt
#
# Additional references:
# http://genomewiki.ucsc.edu/index.php/Chains_Nets
# https://genome.ucsc.edu/goldenPath/help/net.html
# http://genome.ucsc.edu/goldenpath/help/chain.html
# http://asia.ensembl.org/info/website/upload/psl.html
#
# The pipeline depends on GNU parallel, genometools (> 1.5.5) and the following
# tools from UCSC-Kent tookit: faSplit, faToTwoBit, twoBitInfo, blat, axtChain,
# chainSort, chainMergeSort, chainSplit, chainNet, netChainSubset, and liftOver.

require 'yaml'
require 'tempfile'

def add_to_PATH(path)
  return unless path
  return unless File.directory? path
  return if ENV['PATH'].split(':').include? path
  ENV['PATH'] = "#{path}:#{ENV['PATH']}"
end

def to_2bit(fas)
  sh "faToTwoBit #{fas} #{fas.ext('2bit')}"
end

def to_sizes(twobit)
  sh "twoBitInfo #{twobit} stdout | sort -k2nr > #{twobit.ext('sizes')}"
end

def to_ooc(fas)
  sh "blat #{fas} /dev/null /dev/null" \
     " -tileSize=11 -makeOoc=#{run_dir}/#{fas.ext('11.ooc')} -repMatch=100"
end

def extract_cdna(fas, gff)
  sh "gt extractfeat -type exon -join -retainids -coords"                      \
     " -seqfile #{fas} -matchdescstart"                                        \
     " #{gff} > #{gff.ext('.cdna.fa')}"
end

def extract_cds(fas, gff)
  sh "gt extractfeat -type CDS -join -retainids -coords"                       \
     " -seqfile #{fas} -matchdescstart"                                        \
     " #{gff} > #{gff.ext('.cds.fa')}"
end

def extract_pep(fas, gff)
  sh "gt extractfeat -type CDS -translate -join -retainids -coords"            \
     " -seqfile #{fas} -matchdescstart"                                        \
     " #{gff} > #{gff.ext('.pep.fa')}"
end

# Addresses the following issues in the given GFF file.
#
# a. mRNA with subfeatures on different scaffolds. Such annotations are
#    removed.
# b. CDS with no mRNA. An mRNA is added around these CDS.
# c. mRNA with no CDS. Such mRNAs are removed.
def process_gff(gff, out)
  require 'bio/db/gff'

  # Read the lifted gff file into memory and parse it.
  gff3 = Bio::GFF::GFF3.new File.read(gff)

  # Obtain transcripts and their children. genes and other features are not
  # processed.
  transcripts = Hash.new { |h, k| h[k] = [] }
  gff3.records.each do |record|
    # GFF file includes features, comments and directives. We are only
    # interested in "features".
    next unless record.respond_to?(:feature_type)

    # If the feature is a transcript, we consider its ID attribute.
    if record.feature_type =~ /mRNA|transcript/
      transcripts[record.attributes.assoc('ID')] << record
    end

    # If the feature is exon or CDS, we consider its Parent attribute.
    if record.feature_type =~ /exon|CDS/
      transcripts[record.attributes.assoc('Parent')] << record
    end
  end

  transcripts = transcripts.map do |id, annots|
    next unless annots.map(&:seqname).uniq.length == 1    # mRNA on different scaffolds/contigs
    next unless annots.map(&:feature_type).include? 'CDS' # mRNA with no CDS

    # If there's a group of CDS without parent mRNA.
    if !annots.map(&:feature_type).include?('mRNA')
      mrna = [
        annots.first.seqname,
        annots.first.source,
        'mRNA',
        annots.map(&:start).min,
        annots.map(&:end).max,
        nil,
        annots.first.strand,
        nil,
        [["ID", id]]
      ]
      mrna = Bio::GFF::GFF3::Record.new(*mrna)
      annots.unshift(mrna)
    end

    annots
  end.flatten.compact

  tmp = Tempfile.open('lifted')
  tmp.write transcripts.join
  tmp.close
  sh "gt gff3 -tidy -sort -addids -retainids #{tmp.path} > #{out}"
end

def num_sequences(fas)
  `grep '>' #{fas} | wc -l`.strip
end

def num_exact(fas1, fas2)
  Dir.mktmpdir do |dir|
    system "grep -v '>' #{fas1} | sort > #{dir}/#{File.basename fas1}"
    system "grep -v '>' #{fas2} | sort > #{dir}/#{File.basename fas2}"
    comm =
      "comm -12"                                                               \
      " #{dir}/#{File.basename fas1}"                                          \
      " #{dir}/#{File.basename fas2}"                                          \
      " | wc -l"
    `#{comm}`.strip
  end
end

def summarize(source, lifted, outdir)
  File.open("#{outdir}/summary.txt", 'w') do |file|
    %w(cdna cds pep).each do |tag|
      fas1 = source.ext(".#{tag}.fa")
      fas2 = lifted.ext(".#{tag}.fa")
      next unless File.exist?(fas1) || File.exist?(fas2)

      file.puts tag.upcase
      file.puts "  source: #{num_sequences(fas1)}"
      file.puts "  lifted: #{num_sequences(fas2)}"
      file.puts "  exact:  #{num_exact(fas1, fas2)}"
    end
  end
end

def parallel(files, template)
  run_dir = CONFIG[:run_dir]
  name = template.split.first
  jobs = files.map { |file| template % { :this => file } }
  joblst = "#{run_dir}/joblst.#{name}"
  joblog = "#{run_dir}/joblog.#{name}"
  File.write(joblst, jobs.join("\n"))
  sh "parallel --joblog #{joblog} -j #{jobs.length} -a #{joblst}"
end

################################################################################


file "create.liftover" do

  run_dir = CONFIG[:run_dir]
  mkdir "#{run_dir}"

  # keep a copy of setup files
  cp 'Rakefile', "#{run_dir}/"
  cp 'opts.yaml', "#{run_dir}/"

  processes = CONFIG[:processes]
  blat_opts = CONFIG[:blat_opts]

  cp CONFIG[:source_fa], "#{run_dir}/source.fa"
  cp CONFIG[:target_fa], "#{run_dir}/target.fa"

  to_2bit "#{run_dir}/source.fa"
  to_2bit "#{run_dir}/target.fa"

  to_sizes "#{run_dir}/source.2bit"
  to_sizes "#{run_dir}/target.2bit"

  to_ooc "#{run_dir}/source.fa"
  
  # Partition target assembly.
  sh "faSplit sequence #{run_dir}/target.fa #{processes} #{run_dir}/chunk_"

  parallel Dir["#{run_dir}/chunk_*.fa"],
    'faSplit -oneFile size %{this} 5000 %{this}.5k -lift=%{this}.lft &&'       \
    'mv %{this}.5k.fa %{this}'

  # BLAT each chunk of the target assembly to the source assembly.
  parallel Dir["#{run_dir}/chunk_*.fa"],
    "blat -noHead -ooc=#{run_dir}/source.11.ooc #{blat_opts} #{run_dir}/source.fa %{this} %{this}.psl"

  parallel Dir["#{run_dir}/chunk_*.fa"],
    "liftUp -type=.psl -pslQ -nohead"                                          \
    " %{this}.psl.lifted %{this}.lft warn %{this}.psl"

  # Derive a chain file each from BLAT's .psl output files.
  parallel Dir["#{run_dir}/chunk_*.psl.lifted"],
    "axtChain -psl -linearGap=medium"                                          \
    " %{this} #{run_dir}/source.2bit #{run_dir}/target.2bit %{this}.chn"

  # Sort the chain files.
  parallel Dir["#{run_dir}/chunk_*.chn"],
    'chainSort %{this} %{this}.sorted'

  # Combine sorted chain files into a single sorted chain file.
  sh 'chainMergeSort #{run_dir}/*.chn.sorted | chainSplit #{run_dir} stdin -lump=1'
  mv "#{run_dir}/000.chain", "#{run_dir}/combined.chn.sorted"

  # Derive net file from combined, sorted chain file.
  sh "chainNet"                                                                \
     " #{run_dir}/combined.chn.sorted "                                        \
     " #{run_dir}/source.sizes "                                               \
     " #{run_dir}/target.sizes "                                                \
     " #{run_dir}/combined.chn.sorted.net /dev/null"

  # Subset combined, sorted chain file.
  sh "netChainSubset"                                                          \
     " #{run_dir}/combined.chn.sorted.net "                                    \
     " #{run_dir}/combined.chn.sorted'                                         \
     " #{run_dir}/liftover.chn"
end

task 'default' do
  FileUtils.cd Rake.application.original_dir
  fail unless File.exist? 'opts.yaml'

  CONFIG = YAML.load_file 'opts.yaml'
  Array(CONFIG[:add_to_path]).each do |path|
    add_to_PATH path
  end

  # read user provided destination folder
  run_dir = CONFIG[:run_dir]

  Rake.application['create.liftover'].invoke
  
  Array(CONFIG[:lift]).each do |inp|
    outdir = 
      "#{run_dir}/#{File.basename(inp, '.gff3')}-liftover-"                    \
      "#{File.basename(CONFIG[:target_fa], '.fa')}"

    mkdir "#{outdir}"

    out = "#{outdir}/#{outdir}.gff3"

    # Lift over the annotations from source assembly to target assembly.
    sh "liftOver -gff #{inp} #{run_dir}/liftover.chn"                          \
       " #{outdir}/lifted.gff3 #{outdir}/unlifted.gff3"

    process_gff "#{outdir}/lifted.gff3", out

    extract_cdna("#{run_dir}/target.fa", out) if File.exist? inp.ext('cdna.fa')
    extract_cds("#{run_dir}/target.fa", out) if File.exist? inp.ext('cds.fa')
    extract_pep("#{run_dir}/target.fa", out) if File.exist? inp.ext('pep.fa')

    summarize inp, out, outdir
  end
end
