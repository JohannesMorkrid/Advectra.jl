# ------------------------------------------------------------------------------------------
#                                 SMTPClient Extension Tests
# ------------------------------------------------------------------------------------------

# send_mail itself is not called, since the extension loads the user's credentials from
# ext/.env (if present) and would send a real mail

using Test
using Advectra

# Loading warns when there is no ext/.env file, which is expected here
Base.CoreLogging.with_logger(Base.CoreLogging.NullLogger()) do
    @eval using SMTPClient
end
const ext = Base.get_extension(Advectra, :AdvectraSMTPClientExt)

@testset "SMTPClient extension" begin
    @test !isnothing(ext)
    @test hasmethod(send_mail, Tuple{String})

    # Restore the environment afterwards, the .env file sets these variables
    withenv("MAIL_USERNAME" => nothing, "MAIL_PASSWORD" => nothing,
            "MAIL_RECIPIANT" => nothing) do
        mktempdir() do dir
            ext.create_env_file(dir, "me@test.com", "secret", "<you@test.com>")
            @test isfile(joinpath(dir, ".env"))

            ext.load_env_file!(dir)
            @test ENV["MAIL_USERNAME"] == "me@test.com"
            @test ENV["MAIL_PASSWORD"] == "secret"
            @test ENV["MAIL_RECIPIANT"] == "<you@test.com>"

            # Malformed lines and missing files are skipped with a warning
            write(joinpath(dir, ".env"), "# comment\n\nMALFORMED LINE\n")
            @test_logs (:warn, r"Skipping malformed line") ext.load_env_file!(dir)
            @test_logs (:warn, r"No .env file found") ext.load_env_file!(joinpath(dir, "x"))
        end
    end
end
