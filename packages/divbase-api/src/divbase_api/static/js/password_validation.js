/*  
Validates the user password for both registration and password reset: 

1. Checks that the password and confirm password fields match.
2. Checks that the password is not one of the 10,000 most common passwords. 
3. Min password length specified in HTML inputs 

Will prevent form submission if issues. 
*/

const passwordInput = document.getElementById('password');
const confirmPasswordInput = document.getElementById('confirm_password');
const form = document.getElementById("registerForm") || document.getElementById("resetPasswordForm");

let commonPasswords = new Set();

function validatePasswords() {
    const userPassword = passwordInput.value;
    if (commonPasswords.has(userPassword.toLowerCase())) {
        passwordInput.setCustomValidity("That password is one of the 10,000 most common passwords. Please choose a better one.");
    } else {
        // setting validity to empty string means no error. 
        passwordInput.setCustomValidity("");
    }
    if (confirmPasswordInput.value !== userPassword) {
        confirmPasswordInput.setCustomValidity('Passwords do not match');
    } else {
        confirmPasswordInput.setCustomValidity('');
    }
}

async function setupPasswordValidation() {
    try {
        const response = await fetch('/static/files/pwdb_top-10000.json');
        const data = await response.json();
        commonPasswords = new Set(data);

        passwordInput.addEventListener('input', validatePasswords);
        confirmPasswordInput.addEventListener('input', validatePasswords);

        form.addEventListener('submit', (event) => {
            validatePasswords();
            if (!form.checkValidity()) {
                event.preventDefault(); // Stop form submission
            }
        });

    } catch (err) {
        console.error("Failed to load list of common passwords", err);
    }
}

setupPasswordValidation();