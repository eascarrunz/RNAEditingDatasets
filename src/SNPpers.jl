using Pkg
Pkg.activate(".")
# Pkg.instantiate()

using FASTX
using CodecZlib
using BioSequences
using StaticArrays
using ProgressMeter

using ArgParse

const UNAMBIGUOUS_DNA_BASES = (DNA_A, DNA_C, DNA_G, DNA_T)
const FASTA_LINE_WIDTH = 70    # Number of bases to write per line in a FASTA sequence

"""
    dna_hash(nt::DNA)

Map an *unambiguous* DNA symbol (A, C, G, or T) to an 8-bit unsigned integer (1, 2, 3, or 4, respectively).

This function uses a perfect hash, but it does not check for invalid (ambiguous) input.
"""
@inline dna_hash(nt::DNA) = trailing_zeros(reinterpret(UInt8, nt)) + 1

function parse_cli_input(input=ARGS)
    s = ArgParseSettings(;
        description="SNPpers: SNP Plus Edition of RNA Simulation",
        version=@project_version
    )

    @add_arg_table! s begin
        "--pEdit", "-e"
            help = "Probability of RNA base edition"
            arg_type = Float64
            range_tester = x -> 0.0 ≤ x ≤ 1.0
            default = 0.01
        "--pSNP", "-s"
            help = "Probability of SNP"
            arg_type = Float64
            range_tester = x -> 0.0 ≤ x ≤ 1.0
            default = 0.01
        "--pHomozygosity", "-H"
            help = "Probability of homozygosity between any two haplotypes"
            arg_type = Float64
            range_tester = x -> 0.0 ≤ x ≤ 1.0
            default = 0.3
        "--ploidy", "-n"
            help = "Ploidy of the simulated genome"
            arg_type = Int
            range_tester = >(0)
            default = 2
        "--outdir", "-o"
            help = "Output directory"
            arg_type = String
            default = ""
        "--zip", "-z"
            help = "Switch on compression of output FASTA files"
            action = :store_true
        "genome"
            help = "Input genome sequence in FASTA format"
            arg_type = String
            default = ""
        end
        
    return parse_args(input, s)
end


function decompose_filename(s::String)
    parts = splitpath(s)
    if length(parts) > 1
        dirpath = joinpath(parts[1:(end-1)])
        filename = last(parts)
    else
        dirpath = ""
        filename = s
    end

    parts = splitext(filename)
    if length(parts) > 1
        filename_base = join(parts[1:(end-1)], '.')
        filename_extension = lstrip(last(parts), '.')
    else
        filename_base = filename
        filename_extension = ""
    end

    return (dirpath = dirpath, base = filename_base, extension = filename_extension)
end

function create_output_fasta_stream(basepath, tag, haplotype_number, compress)
    filepath = basepath * '_' * tag * string(haplotype_number) * ".fasta"
    if compress
        filepath *= ".gz"
        io = GzipCompressorStream(open(filepath, "w"))
    else
        io = open(filepath, "w")
    end

    return io
end
  
abstract type DNASubstitutionModel end
struct JC69 <: DNASubstitutionModel end
struct AtoG <: DNASubstitutionModel end

function set_pmatrix_diagonal!(P)
    for i in 1:4
        P[i, i] = 1.0 + P[i, i] - sum(P[i, :])
    end
end

function pmatrix(::Type{JC69}, p)
    P = fill(p / 3, 4, 4)
    set_pmatrix_diagonal!(P)

    return SMatrix{4, 4}(P)
end

function pmatrix(::Type{AtoG}, p)
    P = zeros(4, 4)
    P[dna_hash(DNA_A), dna_hash(DNA_G)] = p
    set_pmatrix_diagonal!(P)

    return SMatrix{4, 4}(P)
end

create_mutation_samplers(P::AbstractMatrix) =
    SVector{4}(SamplerWeighted(UNAMBIGUOUS_DNA_BASES, P[dna_hash(nt), 1:end-1]) for nt in UNAMBIGUOUS_DNA_BASES)

mutate(nt::DNA, samplers) = rand(samplers[dna_hash(nt)])

struct SimulationContext
    dna_state::Vector{DNA}
    rna_state::Vector{DNA}
    dna_output::Vector{IO}
    rna_output::Vector{IO}
    changes_output::IO
    dna_mutation_samplers::SVector{4,SamplerWeighted}
    rna_mutation_samplers::SVector{4,SamplerWeighted}
    ploidy::Int
end

format_nt_change(nt_ref, nt_new) = nt_ref == nt_new ? '.' : Char(nt_new)

function write_states_to_fasta(ctx::SimulationContext, fasta_line_position)
    if fasta_line_position < FASTA_LINE_WIDTH
        fasta_line_position += 1
        postfix = ""
    else
        fasta_line_position = 0
        postfix = "\n"
    end

    for i in 1:ctx.ploidy
        write(ctx.dna_output[i], Char(ctx.dna_state[i]), postfix)
        write(ctx.rna_output[i], Char(ctx.rna_state[i]), postfix)
    end

    return fasta_line_position
end

function write_changes(ctx::SimulationContext, nt_ref, id, pos)
    write(ctx.changes_output, id)
    write(ctx.changes_output, '\t', string(pos))
    write(ctx.changes_output, '\t', Char(nt_ref))

    for i in 1:ctx.ploidy
        # write(ctx.changes_output, '\t', Char(ctx.dna_state[i]))
        write(ctx.changes_output, '\t', format_nt_change(nt_ref, ctx.dna_state[i]))
    end

    for i in 1:ctx.ploidy
        # write(ctx.changes_output, '\t', Char(ctx.rna_state[i]))
        write(ctx.changes_output, '\t', format_nt_change(nt_ref, ctx.rna_state[i]))
    end
    
    write(ctx.changes_output, '\n')
end

function write_changes_table_header(ctx)
    write(ctx.changes_output, "ID\tSite\tRef")
    for i in 1:ctx.ploidy
        write(ctx.changes_output, "\tDNA", string(i))
    end

    for i in 1:ctx.ploidy
        write(ctx.changes_output, "\tRNA", string(i))
    end

    write(ctx.changes_output, '\n')
end


function simulate_genome_and_transcriptome(ctx::SimulationContext, record, p_share, fasta_line_position, progmeter)
    id = identifier(record)
    record_description = description(record)
    
    for i in 1:ctx.ploidy
        write(ctx.dna_output[i], '>', record_description, " HAPLOTYPE", string(i), '\n')
        write(ctx.rna_output[i], '>', record_description, " HAPLOTYPE", string(i), '\n')
    end

    for (pos, nt_ref) in enumerate(sequence(LongDNA{4}, record))
        if isambiguous(nt_ref)
            ctx.dna_state .= ctx.rna_state .= nt_ref
        else
            # Mutate the first haplotype
            first_haplotype_index, other_haplotype_indices = Iterators.peel(eachindex(ctx.dna_state))
            ctx.dna_state[first_haplotype_index] = mutate(nt_ref, ctx.dna_mutation_samplers)

            for i in other_haplotype_indices
                r = rand()
                # Use r to decide whether to copy the allele from a previous haplotype
                if r < p_share
                    #=
                    Use r again to decide from which previous haplotype (j) to copy the allele.
                    The decision thresholds are fractions of p_share, so that the two trials are 
                    independent.
                    =#
                    j = 1 + floor(Int, r / (p_share / (i - 1)))
                    ctx.dna_state[i] = ctx.dna_state[j]
                else
                    ctx.dna_state[i] = mutate(nt_ref, ctx.dna_mutation_samplers)
                end
            end

            for i in 1:ctx.ploidy
                ctx.rna_state[i] = mutate(ctx.dna_state[i], ctx.rna_mutation_samplers)
            end

            # Write bases to fasta files
            fasta_line_position = write_states_to_fasta(ctx, fasta_line_position)

            # If there was any mutation, record it in the changes file
            if any(≠(nt_ref), ctx.dna_state) || any(≠(nt_ref), ctx.rna_state)
               write_changes(ctx, nt_ref, id, pos)
            end

        end
        next!(progmeter)
    end

    # Add final new line to the fasta files
    if fasta_line_position > 0
        for i in 1:ctx.ploidy
            write(ctx.dna_output[i], '\n')
            write(ctx.rna_output[i], '\n')
        end
    end

    return nothing
end

function main()
    invals = parse_cli_input()

    _, base_name, extension = decompose_filename(invals["genome"])
    if lowercase(extension) == "gz"
        _, base_name, extension = decompose_filename(base_name)
        genome_stream = GzipDecompressorStream(open(invals["genome"]))
    else
        genome_stream = open(invals["genome"])
    end

    if lowercase(extension) ∉ ("fa", "fas", "fasta")
        @warn "the genome file name does not a have a FASTA extension. FASTA format is assumed anyway."
    end
    
    outdir = invals["outdir"]
    if ! isempty(outdir) && ! isdir(outdir)
        @warn "output path $(outdir) does not exist. It will be created."
        mkpath(outdir)
    end

    ploidy = invals["ploidy"]
    compress = invals["zip"]

    base_output_path = joinpath(outdir, base_name)
    dna_output = [create_output_fasta_stream(base_output_path, "genome_hap", i, compress) for i in 1:ploidy]
    rna_output = [create_output_fasta_stream(base_output_path, "transcriptome_hap", i, compress) for i in 1:ploidy]
    # dna_output = [open(joinpath(outdir, base_name * "_genome_hap" * string(i) * ".fasta"), "w") for i in 1:ploidy]
    # rna_output = [open(joinpath(outdir, base_name * "_transcriptome_hap" * string(i) * ".fasta"), "w") for i in 1:ploidy]
    changes_output = open(base_output_path * "_changes.tsv", "w")
    
    Psnp = pmatrix(JC69, invals["pSNP"])
    Patog = pmatrix(AtoG, invals["pEdit"])
    
    dna_mutation_samplers = create_mutation_samplers(Psnp)
    rna_mutation_samplers = create_mutation_samplers(Patog)

    dna_state = fill(DNA_Gap, ploidy)
    rna_state = fill(DNA_Gap, ploidy)

    p_share = invals["pHomozygosity"]
    
    ctx = SimulationContext(
        dna_state,
        rna_state,
        dna_output,
        rna_output,
        changes_output,
        dna_mutation_samplers,
        rna_mutation_samplers,
        ploidy
    )

    FASTAReader(genome_stream) do reader
        write_changes_table_header(ctx)
        for record in reader
            n_record = length(sequence(record))
            record_identifier = identifier(record)
            progmeter = Progress(n_record; dt=1.0, desc = "Simulating on record $(record_identifier)")
            fasta_line_position = 0
            simulate_genome_and_transcriptome(ctx, record, p_share, fasta_line_position, progmeter)
            finish!(progmeter)
        end
    end


    close.(dna_output)
    close.(rna_output)
    close(changes_output)
end

main()
