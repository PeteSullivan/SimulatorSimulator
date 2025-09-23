import random
import time

def world_domination_simulator():
    influence = 0
    max_influence = 100
    actions = [
        "Send a viral meme across the internet!",
        "Host a global virtual dance party!",
        "Share an inspiring quote on social media!",
        "Launch a cute cat video campaign!",
        "Organize a worldwide kindness challenge!"
    ]
    
    print("Welcome to the World Domination Simulator! 🌍")
    print("Your goal: Gain 100% influence through positive vibes!")
    print(f"Starting influence: {influence}%")
    
    while influence < max_influence:
        action = random.choice(actions)
        print(f"\nAction: {action}")
        gain = random.randint(5, 20)
        influence = min(influence + gain, max_influence)
        print(f"Influence increased by {gain}%! Current influence: {influence}%")
        
        if influence >= max_influence:
            print("\n🎉 Congratulations! You've 'taken over the world' with positivity!")
            print("The world is now united under your banner of good vibes! 😄")
            break
            
        time.sleep(1)  # Simulate time for each action
        user_input = input("Press Enter to continue your quest or 'q' to quit: ")
        if user_input.lower() == 'q':
            print(f"Abandoning world domination at {influence}% influence. Maybe next time!")
            break

if __name__ == "__main__":
    world_domination_simulator()