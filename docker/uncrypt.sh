
function uncrypt_file() {
    local file_to_uncrypt=$1
    local output_file=${2:-${file_to_uncrypt%.*}}
    local password_file=${3:-passwd}
    # shellcheck disable=SC2155
    local gpg_password=$(cat "$password_file")
    gpg --batch --yes --output $output_file --passphrase "${gpg_password}" -d "${file_to_uncrypt}"
}

uncrypt_file $@