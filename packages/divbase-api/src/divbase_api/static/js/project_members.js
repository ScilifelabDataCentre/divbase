// Adds a confirmation dialog for removing a project member
function confirmRemoveMember(userId, userName, projectId) {
    if (confirm(`Are you sure you want to remove ${userName} from this project?`)) {
        const form = document.getElementById('removeMemberForm');
        form.action = `/projects/${projectId}/members/${userId}/remove`;
        form.submit();
    }
}