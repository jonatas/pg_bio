DO $$
BEGIN
    IF NOT EXISTS (SELECT 1 FROM pg_type WHERE typname = 'spectrum') THEN
        CREATE TYPE public.spectrum;
        CREATE OR REPLACE FUNCTION public.spectrum_in(input cstring) RETURNS spectrum LANGUAGE c IMMUTABLE PARALLEL SAFE AS 'pg_bio', 'spectrum_in_wrapper';
        CREATE OR REPLACE FUNCTION public.spectrum_out(input spectrum) RETURNS cstring LANGUAGE c IMMUTABLE PARALLEL SAFE STRICT AS 'pg_bio', 'spectrum_out_wrapper';
        CREATE TYPE public.spectrum (
            INTERNALLENGTH = variable,
            INPUT = public.spectrum_in,
            OUTPUT = public.spectrum_out,
            STORAGE = extended
        );
    END IF;

    IF NOT EXISTS (SELECT 1 FROM pg_type WHERE typname = 'theoreticalfragment') THEN
        CREATE TYPE public.theoreticalfragment;
        CREATE OR REPLACE FUNCTION public.theoreticalfragment_in(input cstring) RETURNS theoreticalfragment LANGUAGE c IMMUTABLE PARALLEL SAFE AS 'pg_bio', 'theoreticalfragment_in_wrapper';
        CREATE OR REPLACE FUNCTION public.theoreticalfragment_out(input theoreticalfragment) RETURNS cstring LANGUAGE c IMMUTABLE PARALLEL SAFE STRICT AS 'pg_bio', 'theoreticalfragment_out_wrapper';
        CREATE TYPE public.theoreticalfragment (
            INTERNALLENGTH = variable,
            INPUT = public.theoreticalfragment_in,
            OUTPUT = public.theoreticalfragment_out,
            STORAGE = extended
        );
    END IF;
END
$$;

CREATE OR REPLACE FUNCTION public.gravy_score(sequence text) RETURNS double precision LANGUAGE c IMMUTABLE PARALLEL SAFE STRICT AS 'pg_bio', 'gravy_score_wrapper';
CREATE OR REPLACE FUNCTION public.isoelectric_point(sequence text) RETURNS double precision LANGUAGE c IMMUTABLE PARALLEL SAFE STRICT AS 'pg_bio', 'isoelectric_point_wrapper';
CREATE OR REPLACE FUNCTION public.create_spectrum(m_z double precision[], intensity double precision[]) RETURNS spectrum LANGUAGE c IMMUTABLE PARALLEL SAFE STRICT AS 'pg_bio', 'create_spectrum_wrapper';
CREATE OR REPLACE FUNCTION public.ppm_difference(mass1 double precision, mass2 double precision) RETURNS double precision LANGUAGE c IMMUTABLE PARALLEL SAFE STRICT AS 'pg_bio', 'ppm_difference_wrapper';
CREATE OR REPLACE FUNCTION public.is_mass_within_ppm(mass1 double precision, mass2 double precision, ppm_tolerance double precision) RETURNS boolean LANGUAGE c IMMUTABLE PARALLEL SAFE STRICT AS 'pg_bio', 'is_mass_within_ppm_wrapper';
CREATE OR REPLACE FUNCTION public.generate_theoretical_fragments(sequence text, max_charge integer, method text) RETURNS TABLE(ion_type text, mass double precision, charge integer) LANGUAGE c IMMUTABLE PARALLEL SAFE STRICT AS 'pg_bio', 'generate_theoretical_fragments_wrapper';
CREATE OR REPLACE FUNCTION public.isotope_fit_score(observed spectrum, theoretical_mz double precision, charge integer, ppm_tolerance double precision) RETURNS double precision LANGUAGE c IMMUTABLE PARALLEL SAFE STRICT AS 'pg_bio', 'isotope_fit_score_wrapper';
CREATE OR REPLACE FUNCTION public.generate_decoy_sequence(sequence text, method text) RETURNS text LANGUAGE c IMMUTABLE PARALLEL SAFE STRICT AS 'pg_bio', 'generate_decoy_sequence_wrapper';
CREATE OR REPLACE FUNCTION public.parse_mgf(filepath text) RETURNS TABLE(title text, pepmass double precision, charge integer, spectrum spectrum) LANGUAGE c STRICT AS 'pg_bio', 'parse_mgf_wrapper';
