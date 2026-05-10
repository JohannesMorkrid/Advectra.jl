module AdvectraSMTPClientExt

using SMTPClient

# Effectively overload by specializing the function call
import Advectra: send_mail
function send_mail(subject::AbstractString; attachment="")
    # Dates.format(now(), RFC1123Format)
    url = "smtps://smtp.gmail.com:465"
    rcpt = split(ENV["MAIL_RECIPIANT"], ",")

    to = ["You"]
    from = "Simulation update"
    message = "Simulation finished" # TODO add better message
    mime_msg = get_mime_msg(message)
    attachments = [attachment]
    if attachment != ""
        body = get_body(to, from, subject, mime_msg; attachments)
    else
        body = get_body(to, from, subject, mime_msg)
    end
    opt = SendOptions(; isSSL=true, username=ENV["MAIL_USERNAME"],
                      passwd=ENV["MAIL_PASSWORD"])
    args = (url, rcpt, from, body, opt) # Have to wrap the args to remove "problem"
    resp = send(args...)
end

"""
    create_env_file(path="", username="mail@provider.com", app_password="password",
    recipiant="<mail@provider.com>"
  Creates an *.env* file in the `path` similar to the `.env.example` with information about 
  the sender: `username` and `app_password` (see https://support.google.com/mail/answer/185833) 
  , and information about the recipiant(s). Pass a string with "\\<mail1>, \\<mail2>, ..." for 
  multiple recipiants.
"""
function create_env_file(path="", username="mail@provider.com", app_password="password",
                         recipiant="<mail@provider.com>")
    filecontent = """
    # Make sure .gitignore works properly as to not share your password!
    MAIL_USERNAME = $username

    # Mail password has to be a app password
    MAIL_PASSWORD = $app_password

    # Mail recipiant
    MAIL_RECIPIANT = $recipiant
    """

    open(joinpath(path, ".env"), "w") do io
        print(io, filecontent)
    end
end

"""
    load_env_file!(file_path=".env")
  Load key-value pairs from a `.env` file in `file_path` into `ENV`.
"""
function load_env_file!(file_path::AbstractString=@__DIR__)
    path = joinpath(file_path, ".env")

    if !isfile(path)
        @warn "No .env file found at $path"
        return nothing
    end

    for line in readlines(path)
        line = strip(line)
        isempty(line) && continue
        startswith(line, "#") && continue

        if occursin("=", line)
            key, value = split(line, "="; limit=2)
            ENV[strip(key)] = strip(value)
        else
            @warn "Skipping malformed line in $path: $line"
        end
    end
end

# Append data into EnvDict
__init__() = load_env_file!()

end