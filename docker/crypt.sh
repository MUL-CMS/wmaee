
function crypt_file() {
    local file_to_crypt=$1
    local output_file=${2:-"${file_to_crypt}.gpg"}
    local password_file=${3:-passwd}
    # shellcheck disable=SC2155
    local gpg_password=$(cat "$password_file")
    gpg --batch --yes --output $output_file --passphrase "${gpg_password}" -c "${file_to_crypt}"
}

# shellcheck disable=SC2068
crypt_file $@